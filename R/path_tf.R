#Functions used for fitting JADE at a series of gammas based on a fit with gamma=0


#Run the full path sequentially
#This is now preferred as of Dec. 2014
#Modified to try to choose better range
jade_path <- function(fit0.file, n.fits, out.file, log.gamma.min=-3, start.step=0.03, tol=1e-5, mean=NULL, max_it=10000, log.gamma.max=20, max.iter=NULL, buffer=0.01, temp.file=NULL, restart.file=NULL){

		if(is.null(temp.file)){
				z <- unlist(strsplit(out.file, ".RData"))[1]
				temp.file <- paste(z, "_temp.RData", sep="")
		}
	#Set up for fits
	fit0 <- getobj(fit0.file)
	log.gammas <- c(-Inf, log.gamma.min)

	p <- dim(fit0$y)[1]
	K <- dim(fit0$y)[2]
	l1.total0 <- 0
	sep.total0 <- 0
	for(j in 1:(K-1)){
    for(l in (j+1):K){
			sep.total0 <- sep.total0 + sum(abs(fit0$fits[,j]-fit0$fits[,l]) > tol)
			l1.total0 <- l1.total0 + sum(abs(fit0$fits[,j] - fit0$fits[,l]))
		}
	}
	cat("K: ", K, " p: ", p, " l1.total0: ", l1.total0, " sep.total0: ", sep.total0, "\n")

	l1.gap <- l1.total0/n.fits
	#cat("Gap: ", l1.gap, "\n")
	if(is.null(max.iter)) max.iter <- n.fits*10
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

	theta_init=fit0$fits
	done <- FALSE

	if(!is.null(restart.file)){
			path.temp <- getobj(restart.file)
			i <- length(path.temp$gammas)
			log.gammas <- log10(path.temp$gammas)
			sep <- path.temp$sep
			l1.total <- path.temp$l1.total
			tol <- path.temp$tol
			fits <- path.temp$JADE_fits
			sep.total <- path.temp$sep.total
			new.gamma <- log.gammas[-1]
			closest.idx <- which.min(abs(log.gammas[-1]-new.gamma))+1
			cat("Theta init idx ", closest.idx, "\n")
			theta_init <- fits[[closest.idx]]$fits
	}
	while(!done){
		g <- 10^log.gammas[i]
		cat("Gamma", g, "\n")
		fit <- jade_multi_tf_admm(y=fit0$y, gamma=g, lambda=fit0$lambda, sample.size=fit0$sample.size, sds=fit0$sds, ord=fit0$ord, positions=fit0$pos, scale.pos=fit0$scale.pos, theta_init=theta_init, fit_var=fit0$fit_var, max_it=max_it, tol=tol)
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
		cat(i, "log.gamma: ", log.gammas[i], "n rep: ", fit$n, " converged: ", fit$converged, " sep.total: ", sep.total[i], " l1.total: ", l1.total[i], "\n")
		converged[i] <- fit$converged


		#Update the gap size
		l1.top <- min(l1.total[sep.total >= sep.total0])
		l1.gap <- l1.top/n.fits
		lg.top <- max(log.gammas[sep.total >= sep.total0])
		#l1.maxgap <- 2*l1.gap


		#Try to find the next gamma to evaluate
		if(i >=6){
			#j <- order(c(l1.total[l1.total <= l1.top], 0), decreasing=TRUE)
			#lt <- c(l1.total[l1.total <= l1.top], 0)[j]
			lt <- sort(c(l1.total[ l1.total <= l1.top], 0), decreasing=TRUE)

			#If we are still at the top somewhere
			if(l1.top == min(l1.total)){
				cat("Top!\n")
				l1.target <- min(l1.total) - l1.gap
				new.gamma <- find_new_gamma_tf(l1.target=l1.target, lg.top=lg.top, l1.total=l1.total[-1], log.gammas=log.gammas[-1], l1.gap=l1.gap, buffer=buffer)
				if(new.gamma$warn){
					cat("Hmmm. Check what is happening\n")
					new.gamma$new.gamma <- max(log.gammas) + buffer
				}
				new.gamma <- new.gamma$new.gamma
			}else if(all(-1*diff(lt) <= (l1.gap*2)) | i >= max.iter){
				cat("Bottom?\n")
				#If the path is dense enough
				#And we got close to the bottom
				if( min(l1.total) < (tol*p) | min(sep.total)==0 | max(log.gammas) ==log.gamma.max | i >= max.iter){
					done <- TRUE
					cat("finishing\n")
					if(i==max.iter) cat("Warning: Path may not be complete!\n")
					break
				}else{
					#If the path is dense but we didn't get close to the bottom
					l1.target <- (min(l1.total)-0)/2
					new.gamma <- find_new_gamma_tf(l1.target=l1.target, lg.top=lg.top, l1.total=l1.total[-1], log.gammas=log.gammas[-1], l1.gap=l1.gap, buffer=buffer)
					#if(is.na(new.gamma)) new.gamma <- max(log.gammas)+1
					new.gamma <- new.gamma$new.gamma
				}
			}else{
				cat("Middle\n")
				#If we are in the middle but the path isn't dense enough
				hole.idx <- which(-1*diff(lt) > l1.gap)
				l1.target <- lt[hole.idx]-l1.gap
				new.gamma <- find_new_gamma_tf(l1.target=l1.target, lg.top=lg.top, l1.total=l1.total[-1], log.gammas=log.gammas[-1], l1.gap=l1.gap, buffer=buffer)
				if(new.gamma$warn & (min(l1.total) < (tol*p) | min(sep.total) == 0)){
						done <- TRUE
						cat("finishing\n")
						cat("Warning: Path may have holes!\n")
						break
				}
				new.gamma <- new.gamma$new.gamma
			}
		}else{
			#if i==2,3,4
			new.gamma <- log.gammas[i] + start.step
			if((l1.total0 - l1.total[i])  > l1.gap | sep.total[i] < sep.total0){
				cat("Warning: May desire to start path earlier\n")
				sep.total0 <- sep.total[i]
			}
		}
		closest.idx <- which.min(abs(log.gammas[-1]-new.gamma))+1
		cat("Theta init idx ", closest.idx, "\n")
		theta_init <- fits[[closest.idx]]$fits

		log.gammas <- c(log.gammas, min(new.gamma, log.gamma.max))
		i <- i+1
		cat("Next: ", log.gammas[i], " ")
		path.temp <- list("JADE_fits"=fits, "sep"=sep, "l1.total" = l1.total, "sep.total"=sep.total, "gammas"=10^(log.gammas), "tol"=tol)
		save(path.temp, file=temp.file)
	}

	bic <- lapply(fits[-1], FUN=jade_bic, tol_beta=tol)
	bic <- matrix(unlist(bic), nrow=length(fits)-1, byrow=TRUE)
	path <- list("JADE_fits"=fits, "sep"=sep, "l1.total" = l1.total, "sep.total"=sep.total, "gammas"=10^(log.gammas), "bic"=bic[,1], "df"=bic[,2], "tol"=tol)

	if(!is.null(mean)){
		path$SSES <- unlist(lapply(path$JADE_fits, FUN=jade_sse, mean=mean))
		path$WT_SSES <- unlist(lapply(path$JADE_fits, FUN=jade_wtd_sse, mean=mean))
	}

	cat("Done!\n")
		unlink(temp.file)
	save(path, file=out.file)
}

find_new_gamma_tf <- function(l1.target, lg.top, l1.total, log.gammas, buffer, l1.gap){

	j <- length(l1.target)

	y=l1.total[order(log.gammas)]
	x=sort(log.gammas)
	if(length(x) < 15) ord <- 1
		else ord <- 2
	if(length(x) < 25) k <- 3
		else k <- 5

	cat(length(x), ord, k, "\n")
	#cat(y, "\n", x, "\n")
	x.new <- seq(-18, 20, by=(buffer/2))
	sm <- genlasso:::trendfilter(pos=x, y=y, ord=ord)
	sm.cv <- cv.trendfilter(sm, k=k)
	z <- coef.genlasso(sm, lambda=sm.cv$lambda.1s)
	zz = .Call("tf_predict_R",
                    sBeta = as.double(z$beta),
                    sX = as.double(x),
                    sN = length(y),
                    sK = as.integer(ord),
                    sX0 = as.double(x.new),
                    sN0 = length(x.new),
                    sNLambda = 1,
                    sFamily = 0,
                    sZeroTol = as.double(1e-6),
                    package = "glmgen")
	#lines(x.new, zz)
	new.gamma <- NA
	i <- 1
	warn <- FALSE
	while(is.na(new.gamma) & i <=j){
		new.gamma <- x.new[ which.min(abs(zz-l1.target[i]))]
		cat(new.gamma, " ")
		new.gamma <- sort(c(lg.top-(2*buffer), new.gamma, lg.top+1))[2]
		ng <- new.gamma
		cat(new.gamma, " ")
		#Don't repeat a gamma we have already tried
		while(min(abs(log.gammas-new.gamma)) <= buffer) new.gamma <- new.gamma+(buffer/2)
		cat(new.gamma, "\n")
		if(abs(new.gamma-ng) > 15*buffer){
			if(i < j) new.gamma <- NA
			else warn <- TRUE
		}
		i <- i+1
	}
	return(list("new.gamma"=new.gamma, "warn"=warn))
}
