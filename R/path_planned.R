#Uses estimated start and stop points: log.gamma.start, log.gamma.stop. These should be good estimates
#Rather than checking l1.gap we will just move from start to stop in even steps
#Hard stop means go all the way to log.gamma.stop and don't go past even if path isn't done
#' Fit JADE at a sequence of gamma values.
#' @description This function requires having previously fit the data at \code{gamma}=0. JADE will be fit
#' at a serries of \code{gamma} values evenly spaced on the log scale. The spacing is determined by
#' a starting point, a stopping point and the number of fits desired.
#'
#' The start and end points should be pretty good guesses but don't need to be exact.
#' If \code{log.gamma.start} is too large \code{jade_path_planned} will back up until it finds a value of
#' \code{gamma} with no fusion between profiles. Similarly if \code{log.gamma.stop} is too small and the path isn't complete,
#' the function will keep going until all the profiles are fused.\code{log.gamma.stop}.
#' Guessing the start and stop locations can have a few negative consequences:
#'
#' -If \code{log.gamma.start} is much too large it may take a (prohibitively) long time
#' to fit the initial fit or convergence might not be reached.
#' Warm starts help a lot with convergence times so starting in the middle of the path can be costly.
#'
#' - If \code{log.gamma.start} and \code{log.gamma.stop} are way too far apart, the chosen step size will be too large
#' which will have a similar effect as starting in the middle of the path. It will also result in a sparsely populated
#' path which might not be very useful.
#'
#' - If the endpoints are way too close together the step size will be
#' way small and many more fits than \code{n.fits} will be calculated.
#' @param  fit0 Either a jade object or a file containing a jade object produced by fitting with \code{gamma=0}
#' @param out.file Name of a file to save the results to.
#' @param log.gamma.start Starting value of \code{log10(gamma)}.
#' @param n.fits Desired number of fits along the path. Actual results may vary.
#' \code{n.fits} can only be missing if \code{step.size} is provided.
#' @param log.gamma.stop Speculative stopping point. If \code{log.gamma.stop} and \code{n.fits}
#' are provided \code{step.size = (log.gamma.stop-log.gamma.start)/n.fits}.
#' @param step.size Alternative to specifying \code{log.gamma.stop} and \code{n.fits}
#' @param temp.file Name a temp file. If missing will use a default based on \code{out.file}.
#' This file will be deleted at the end.
#' @param max.it Maximum number of jade iterations.
#' @param log.gamma.max Absolute largest value of \code{gamma} that will be fit.
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
jade_path_planned <- function(fit0, out.file, log.gamma.start,
                              n.fits=NULL, log.gamma.stop=NULL, step.size=NULL,
                              tol=1e-3, verbose = TRUE,
                              temp.file=NULL, max.it=10000, log.gamma.max=20,
                              restart.file=NULL, hard.stop=FALSE){


	if(is.null(temp.file)){
    z <- unlist(strsplit(out.file, ".RData"))[1]
		temp.file <- paste(z, "_temp.RData", sep="")
	}

  if((is.null(log.gamma.stop) | is.null(n.fits) )& is.null(step.size)){
    stop("Specify either log.gamma.stop and n.fits or specify step.size.")
  }else if(is.null(step.size)){
    step.size <- (log.gamma.stop-log.gamma.start)/n.fits
  }

  if(class(fit0)=="character"){
    fit0.file <- fit0
    fit0 <- getobj(fit0.file)
  }
  alg <-  fit0$algorithm
  sep0 <- get_sep(fit0$fits, tol)

  #Set up
	log.gammas <- c(-Inf, log.gamma.start)

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

	if(verbose) cat("K: ", K, " p: ", p, " l1.total0: ", l1.total0, " sep.total0: ", sep.total0, "\n")

	#Structures for path. The first element of fits and l1 and sep lists are for gamma=0
	#sep is a list of lists of  matrices. The first collumn is for gamma=0
	sep <- list()
  for(j in 1:(K-1)){
		sep[[j]] <- list()
    for(l in (j+1):K){
			sep[[j]][[l-j]] <- matrix(sep0[[j]][[l-j]], nrow=p, ncol=1)
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
	  u.beta0=NULL
	}else if(alg=="gd"){
	  duals0=NULL
	}
	done <- FALSE

	#Restarting from a temp file
	if(!is.null(restart.file)){
			path.temp <- getobj(restart.file)
			i <- length(path.temp$gammas)
			log.gammas <- log10(path.temp$gammas)
			sep <- path.temp$sep
			l1.total <- path.temp$l1.total
			tol <- path.temp$tol
			fits <- path.temp$JADE_fits
			sep.total <- path.temp$sep.total
			new.gamma <- log.gammas[i]
			closest.idx <- which.min(abs(log.gammas[-i][-1]-new.gamma))+1
			#cat("Theta init idx ", closest.idx, "\n")
			theta0 <- fits[[closest.idx]]$fits
			if(alg=="admm"){
			  #u.alpha0=fits[[closest.idx]]$u.alpha
			  #u.beta0=fits[[closest.idx]]$u.beta
			}else if(alg=="gd"){
			  duals0=fits[[closest.idx]]$duals
			}
	}

	while(!done){
		g <- 10^log.gammas[i]
		#if(verbose) cat("Gamma", g, "\n")
		if(alg=="admm"){
		  fit <- jade_admm(y=fit0$y, gamma=g, pos=fit0$pos, scale.pos=fit0$scale.pos,
		                   lambda=fit0$lambda, sample.size=fit0$sample.size, ord=fit0$ord,
		                   sds=fit0$sds, fit.var=fit0$fit.var,
		                   theta0=theta0, u.alpha0=u.alpha0, u.beta0=u.beta0,
		                   tol=tol,  max.it=max.it)
		}else if(alg=="gd"){
		  fit <- jade_gd(y=fit0$y, gamma=g, pos=fit0$pos, scale.pos=fit0$scale.pos,
		                 lambda1=fit0$lambda1, lambda2=fit0$lambda2,
		                 sample.size=fit0$sample.size, ord=fit0$ord,
		                 sds=fit0$sds, fit.var=fit0$fit.var,
		                 theta0=theta0, duals0=duals0,
		                 sep.tol=tol,  max.it=max.it)
		}
		fits[[i]] <- fit

		#Add to path elements: separation matrix, converged, l1.total, raw fits
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
		if(verbose) cat(i, "log.gamma: ", log.gammas[i], "n rep: ", fit$n, " converged: ", fit$converged, " sep.total: ", sep.total[i], " l1.total: ", l1.total[i], "\n")
		converged[i] <- fit$converged

		#Find the next gamma to evaluate
		if(!hard.stop){
			if((log.gammas[i] - log.gamma.start) < 2*step.size & sep.total[i] < (0.95*sep.total0)){
			  #Near the top - check if we started too far in:
			  new.gamma <- min(log.gammas[-1])-step.size
			}else if(sep.total[i] == 0 | l1.total[i] < tol | log.gammas[i] >= log.gamma.max){
			  #Time to stop
			  done <- TRUE
		  	new.gamma <- max(log.gammas[-1]) + step.size
			}else{
				new.gamma <- max(log.gammas[-1]) + step.size
			}
		}else{
			new.gamma <- max(log.gammas[-1]) + step.size
			if(log.gammas[i] > log.gamma.stop){
				done <- TRUE
			}
		}
		#Find the fit with the closest gamma to the next value.
		#Use solutions from this fit as new theta0 value.
		closest.idx <- which.min(abs(log.gammas[-1]-new.gamma))+1
		#cat("Theta init idx ", closest.idx, "\n")
		theta0 <- fits[[closest.idx]]$fits
		if(alg=="admm"){
		  #u.alpha0 <- fits[[closest.idx]]$u.alpha
		  #u.beta0 <- fits[[closest.idx]]$u.beta
		}else if(alg=="gd"){
		  duals0 <-  fits[[closest.idx]]$duals
		}
		log.gammas <- c(log.gammas, new.gamma)
		i <- i+1
		path.temp <- list("JADE_fits"=fits, "sep"=sep, "l1.total" = l1.total, "converged"=converged,
		                  "sep.total"=sep.total, "gammas"=10^(log.gammas), "tol"=tol)
		save(path.temp, file=temp.file)
		#cat("Next: ", log.gammas[i], " ")
	}

	#bic <- lapply(fits[-1], FUN=jade_bic, tol_beta=tol)
	#bic <- matrix(unlist(bic), nrow=length(fits)-1, byrow=TRUE)

	log.gammas <- log.gammas[ -(length(log.gammas))]
	path <- list("JADE_fits"=fits, "sep"=sep, "converged"=converged,
	             "l1.total" = l1.total, "sep.total"=sep.total,
	             #"bic"=bic[,1], "df"=bic[,2],
	             "gammas"=10^(log.gammas), "tol"=tol)
	if(verbose) cat("Done!\n")
	save(path, file=out.file)
	unlink(temp.file)
	return(path)
}


