#Functions used for fitting JADE at a series of gammas based on a fit with gamma=0
#Run the full path sequentially
#This is now preferred as of Dec. 2014
#Modified to try to choose better range

#' Fit JADE at a sequence of gamma values without knowing much about where to start.
#' @description This function requires having previously fit the data at \code{gamma}=0. JADE will be fit
#' at a serries of \code{gamma} values The function tries to fit JADE at a series of \code{gamma}
#' values so that fits occur at evenly spaced values of \code{l1.total}.
#'
#' @param fit0 Either a jade object or a file containing a jade object produced by fitting with \code{gamma=0}
#' @param out.file Name of a file to save the results to.
#' @param n.fits Desired number of fits along the path. Actual results may vary.
#' @param temp.file Name a temp file. If missing will use a default based on \code{out.file}.
#' This file will be deleted at the end.
#' @param max.it Maximum number of jade iterations.
#' @param max.fits Maximum number of fits.
#' @param log.gamma.minlog.gamma.max Smallest and largest values of \code{gamma} that are allowed.
#' @param restart.file Provide if restarting from a temporary file.
#' @param hard.stop Stop at \code{log.gamma.stop} regardless of if the path is done or not.
#'
#' @return The output of this function is both returned and saved to a file. Partial results
#' are saved to a temporary file along the way. The path objects is a list with several components:
#' \describe{
#'  \item{\code{gammas}}{A list of length n of gamma values at which JADE was fit. The first element
#'  is always zero.}
#'   \item{\code{JADE_fits}}{A list of n JADE objects in the same order as \code{gammas}.}
#'   \item{\code{l1.total}}{Vector of length n giving the total L1 distance between pairs of profiles at each value of \code{gamma}}
#'   \item{\code{sep.total}}{Vector of length n giving the total number of separated sites between all pairs of profiles}
#'   \item{\code{sep}}{List of lists of matrices giving the pairwise separation
#'   between profiles. sep[[i]][[j-i]] is a \eqn{p \times n} matrix which describes the
#'   separation between profiles for group i and group j.}
#'   \item{\code{tol}}{The tolerance at which the seaparation in \code{sep.total} and \code{sep}
#'   were calculated}
#' }
#' @export
jade_path <- function(fit0, n.fits, out.file, temp.file=NULL,
                            max.it=10000, log.gamma.min=-3, log.gamma.max=20,
                            start.step=0.03, tol=1e-3, max.fits= 10*n.fits,
                            buffer=0.001, restart.file=NULL, verbose=TRUE,
                            adjust.rho.alpha=TRUE){


	if(is.null(temp.file)){
    z <- unlist(strsplit(out.file, ".RData"))[1]
    temp.file <- paste(z, "_temp.RData", sep="")
	}


  if(class(fit0)=="character"){
    fit0.file <- fit0
    fit0 <- getobj(fit0.file)
  }
  alg = fit0$algorithm
  sep0 <- get_sep(fit0$fits, tol)

	#Set up for fits
	log.gammas <- c(-Inf, log.gamma.min)

	p <- dim(fit0$y)[1]
	K <- dim(fit0$y)[2]
	l1.total0 <- 0
	sep.total0 <- 0
	for(j in 1:(K-1)){
    for(l in (j+1):K){
			sep.total0 <- sep.total0 + sum(sep0[[j]][[l-j]])
			l1.total0 <- l1.total0 + sum(abs(fit0$fits[,j] - fit0$fits[,l]))
		}
	}
	cat("K: ", K, " p: ", p, " l1.total0: ", l1.total0, " sep.total0: ", sep.total0, "\n")

	l1.gap <- l1.total0/n.fits
	#cat("Gap: ", l1.gap, "\n")

	#l1.maxgap <- 2*l1.gap

	#Structures for path
	sep <- list()
  for(j in 1:(K-1)){
		sep[[j]] <- list()
    for(l in (j+1):K){
			sep[[j]][[l-j]] <- matrix(NA, nrow=p, ncol=0)
    }
	}
	converged <- c(TRUE)
	l1.total <- c(l1.total0)
	sep.total <- c(sep.total0)
	fits <- list()
	fits[[1]] <- fit0

	i <- 2

	theta0=fit0$fits

	if(alg=="admm"){
	  u.alpha0=NULL
	  rho.alpha=NULL
	  u.beta0=NULL
	  rho.beta=1
	}else if(alg=="gd"){
	  duals0=NULL
	}
	done <- FALSE

	#Restarting from a temp file
	if(!is.null(restart.file)){
			path.temp <- getobj(restart.file)
			i <- length(path.temp$JADE_fits)+1
			log.gammas <- log10(path.temp$gammas)
			sep <- path.temp$sep
			l1.total <- path.temp$l1.total
			tol <- path.temp$tol
			fits <- path.temp$JADE_fits
			sep.total <- path.temp$sep.total
      converged <- path.temp$converged
      if(length(log.gammas)==i){
        new.gamma <- log.gammas[i]
        log.gammas <- log.gammas[-i]
      }else{
			  #Find next gamma value
        if(all(sep.total[-1] < sep.total0)) sep.total0 <- sep.total[2]
			  l1.top <- min(l1.total[sep.total >= sep.total0 & is.finite(log.gammas)])
			  l1.gap <- l1.top/n.fits
			  lg.top <- max(log.gammas[sep.total >= sep.total0])
			  if(!is.finite(lg.top)) lg.top <- min(log.gammas[-1])
			  keep.fits <- which(l1.total <= l1.total0 & is.finite(log.gammas))
			  if(length(keep.fits) < 6) keep.fits <- which(is.finite(log.gammas))
			  new.gamma <- find_new_gamma(l1.total=l1.total[keep.fits], sep.total=sep.total[keep.fits],
			                            log.gammas=log.gammas[keep.fits], start.step=start.step,
			                            l1.gap=l1.gap, l1.top=l1.top,
			                            tol=tol, buffer=buffer, lg.top = lg.top)
      }
			closest.idx <- which.min(abs(log.gammas[-1]-new.gamma))+1
			#cat("Theta init idx ", closest.idx, "\n")
			theta0 <- fits[[closest.idx]]$fits
			if(alg=="admm"){
			  u.alpha0=fits[[closest.idx]]$u.alpha
			  rho.alpha=fits[[closest.idx]]$rho.alpha
			  u.beta0=fits[[closest.idx]]$u.beta
			  rho.beta=fits[[closest.idx]]$rho.beta
			}else if(alg=="gd"){
			  duals0=fits[[closest.idx]]$duals
			}
			log.gammas <- c(log.gammas, min(new.gamma, log.gamma.max))
	}

	while(!done){
		g <- 10^log.gammas[i]
		#cat("Gamma", g, "\n")
		if(alg=="admm"){
		  fit <- jade_admm(y=fit0$y, gamma=g, pos=fit0$pos, scale.pos=fit0$scale.pos,
		                   lambda=fit0$lambda, sample.size=fit0$sample.size, ord=fit0$ord,
		                   sds=fit0$sds, fit.var=fit0$fit.var,
		                   theta0=theta0, u.alpha0=u.alpha0, u.beta0=u.beta0, rho.beta=rho.beta,
		                   rho.alpha=rho.alpha, tol=tol,  max.it=max.it,
		                   adjust.rho.alpha=adjust.rho.alpha)
		}else if(alg=="gd"){
		  fit <- jade_gd(y=fit0$y, gamma=g, pos=fit0$pos, scale.pos=fit0$scale.pos,
		                   lambda1=fit0$lambda1, lambda2=fit0$lambda2,
		                   sample.size=fit0$sample.size, ord=fit0$ord,
		                   sds=fit0$sds, fit.var=fit0$fit.var,
		                   theta0=theta0, duals0=duals0,
		                   sep.tol=tol,  max.it=max.it)
		}

		fits[[i]] <- fit

		#Path elements: separation matrix, converged, l1.total, raw fits
		sep.total <- c(sep.total, 0)
		l1.total <- c(l1.total, 0)
		converged <- c(converged, 0)
  	for(j in 1:(K-1)){
    	for(l in (j+1):K){
				sep[[j]][[l-j]] <- cbind( sep[[j]][[l-j]], fit$sep[[j]][[l-j]])
				sep.total[i] <- sep.total[i] + sum(fit$sep[[j]][[l-j]])
				l1.total[i] <- l1.total[i] + sum(abs(fit$fits[,j] - fit$fits[,l]))
			}
		}
		if(verbose) cat(i, "log.gamma: ", log.gammas[i], "n rep: ", fit$n,
		                " converged: ", fit$converged, " sep.total: ", sep.total[i],
		                " l1.total: ", l1.total[i], "\n")
		converged[i] <- fit$converged


		if(i == 2 & sep.total[i] < sep.total0){
		  cat("Warning: May need to start path earlier")
		  sep.total0 <- sep.total[i]
		}
		#Update the gap size
		#We don't care about the density of the path when sep.total >= sep.total0.
		#l1.top is the smallest value of l1.total achieved while the profiles are
		#still as separated as they are with no penalty.
		#We want a dense path between l1.top and 0
		l1.top <- min(l1.total[sep.total >= sep.total0 & is.finite(log.gammas)])
		l1.gap <- l1.top/n.fits
		lg.top <- max(log.gammas[sep.total >= sep.total0])

    if(!is.finite(lg.top)) lg.top <- min(log.gammas[-1])

		#Find the next gamma to evaluate
		keep.fits <- which(l1.total <= l1.total0 & is.finite(log.gammas))
		if(length(keep.fits) < 6) keep.fits <- which(is.finite(log.gammas))
		new.gamma <- find_new_gamma(l1.total=l1.total[keep.fits], sep.total=sep.total[keep.fits],
		                            log.gammas=log.gammas[keep.fits], start.step=start.step,
		                            l1.gap=l1.gap, l1.top=l1.top,
		                            tol=tol, buffer=buffer, lg.top = lg.top)

		if(is.na(new.gamma) | max(log.gammas) ==log.gamma.max | i >= max.fits){
		  if(verbose) cat("Finishing\n")
		  done <- TRUE
		  break
		}

		#Find the fit with the closest gamma to the next value.
		#Use solutions from this fit as new theta0 value.
		closest.idx <- which.min(abs(log.gammas[-1]-new.gamma))+1
		#cat("Theta init idx ", closest.idx, "\n")
		theta0 <- fits[[closest.idx]]$fits
		if(alg=="admm"){
		  u.alpha0=fits[[closest.idx]]$u.alpha
		  rho.alpha=fits[[closest.idx]]$rho.alpha
		  u.beta0=fits[[closest.idx]]$u.beta
		  rho.beta=fits[[closest.idx]]$rho.beta
		}else if(alg=="gd"){
		  duals0 <-  fits[[closest.idx]]$duals
		}
		log.gammas <- c(log.gammas, min(new.gamma, log.gamma.max))
		i <- i+1
		#cat("Next: ", log.gammas[i], " ")
		path.temp <- list("JADE_fits"=fits, "sep"=sep, "l1.total" = l1.total, "converged"=converged,
		                  "sep.total"=sep.total, "gammas"=10^(log.gammas), "tol"=tol)
		save(path.temp, file=temp.file)
	}

	path <- list("JADE_fits"=fits, "sep"=sep, "l1.total" = l1.total, "converged"=converged,
		             #"bic"=bic[,1], "df"=bic[,2],
		             "sep.total"=sep.total, "gammas"=10^(log.gammas), "tol"=tol)

	cat("Done!\n")
	save(path, file=out.file)
	unlink(temp.file)
  return(path)
}
