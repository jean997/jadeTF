
get_cv_err <- function(orig.path.file, cv.file.list=NULL, n.folds=10, loss.type="mse", use.converged.only=TRUE, control.l1=TRUE){

	stopifnot(loss.type %in% c("mse", "l1"))

	if(is.null(cv.file.list)){
		file.start <- unlist(strsplit(orig.path.file,  split=".RData", fixed=TRUE))[1]
		cv.file.list <- paste(file.start, 1:n.folds, "RData", sep=".")
	}

	#Path with full data
	orig.path <- getobj(orig.path.file)
	orig.y <- orig.path$JADE_fits[[1]]$y
	orig.na <- which(is.na(orig.y))
	K <- dim(orig.y)[2]
	p <- dim(orig.y)[1]
	n.gamma <- length(orig.path$JADE_fits)

	#Correction for previous bug in path_planned which resulted in two coppies of the last fit
	if(length(orig.path$gammas) == n.gamma + 1) orig.path$gammas <- orig.path$gammas[-(n.gamma + 1)]


	#Number of separated sites
	if(is.null(orig.path$sep.total)){
		orig.path$sep.total <- rep(0, n.gamma)
		for(j in 1:(K-1)){
			for(l in (j+1):K){
				orig.path$sep.total <- orig.path$sep.total + colSums(orig.path$sep[[j]][[l-j]])
			}
		}
	}

	#Track l1 distance between curves
	if(is.null(orig.path$l1.total)){
		orig.path$l1.total <- unlist(lapply(orig.path$JADE_fits, FUN=function(x, K){
					z <- 0
					for(j in 1:(K-1)){
						for(l in (j+1):K){
							z <- z+ sum(abs(x$fits[,j]-x$fits[,l]))
						}
					}
					return(z)
				}, K=K))
	}

	#Track convergence
	converged <- list()
	converged[[1]] <- unlist(lapply(orig.path$JADE_fits, FUN=function(fit){
			if(is.null(fit$converged)) return(TRUE)
			return(fit$converged)
			}))

	keep.fits <- list(rep(TRUE, n.gamma))
	if(use.converged.only){
		keep.fits[[1]][!converged[[1]]]<- FALSE
	}
	if(control.l1){
		keep.fits[[1]][ orig.path$l1.total > orig.path$l1.total[1]] <- FALSE
	}
	#CV by matching gamma
	cv.err.gamma <- matrix(0, nrow=n.folds, ncol=sum(keep.fits[[1]]))
	#CV by matching l1 distance
	cv.err.l1 <-  matrix(0, nrow=n.folds, ncol=sum(keep.fits[[1]]))

	log.gamma.list <- list(log10(orig.path$gammas))
	l1.list <- list(orig.path$l1.total)
	sep.list <- list(orig.path$sep.total)

	n.test <- c()

	#Collect info from each fold
	for(i in 1:n.folds){
		path <- getobj(cv.file.list[i])
		cv.y <- path$JADE_fits[[1]]$y

		#Correction for bug in path_planned
		if(length(path$gammas) == length(path$JADE_fits) + 1) path$gammas <- path$gammas[-(length(path$JADE_fits) + 1)]

		#sep.total
		if(is.null(path$sep.total)){
			path$sep.total <- rep(0, length(path$gammas))
			for(j in 1:(K-1)){
				for(l in (j+1):K){
					path$sep.total <- path$sep.total + colSums(path$sep[[j]][[l-j]])
				}
			}
		}

		#l1.total
		if(is.null(path$l1.total)){
			path$l1.total <- unlist(lapply(path$JADE_fits, FUN=function(x, K){
					z <- 0
					for(j in 1:(K-1)){
						for(l in (j+1):K){
							z <- z+ sum(abs(x$fits[,j]-x$fits[,l]))
						}
					}
					return(z)
				}, K=K))
		}

		#convergence
		converged[[(i+1)]] <- unlist(lapply(path$JADE_fits, FUN=function(fit){
			if(is.null(fit$converged)) return(TRUE)
			return(fit$converged)
			}))

		keep.fits[[(i+1)]] <- rep(TRUE, length(path$gammas))
		if(use.converged.only){
			keep.fits[[(i+1)]][!converged[[(i+1)]] ]<- FALSE
		}
		if(control.l1){
			keep.fits[[(i+1)]][ path$l1.total > path$l1.total[1]] <- FALSE
		}

		log.gamma.list[[i+1]] <- log10(path$gammas)
		l1.list[[i+1]] <- path$l1.total
		sep.list[[i+1]] <- path$sep.total

		cv.na <- which(is.na(cv.y))
		test.idx <- cv.na[!cv.na %in% orig.na]
		n.test <- c(n.test, length(test.idx))
		if(loss.type=="mse"){
			cv.err <- unlist( lapply(path$JADE_fits, FUN=function(x, orig.y, test.idx){
					sum((x$fits[test.idx]-orig.y[test.idx])^2)
					}, orig.y=orig.y, test.idx=test.idx))/n.test[i]
		}else if(loss.type=="l1"){
			cv.err <- unlist( lapply(path$JADE_fits, FUN=function(x, orig.y, test.idx){
					sum(abs(x$fits[test.idx]-orig.y[test.idx]))
					}, orig.y=orig.y, test.idx=test.idx))/n.test[i]
		}
		cv.err.l1[i,] <- approx(x=path$l1.total[keep.fits[[(i+1)]]], y=cv.err[keep.fits[[(i+1)]]], xout=orig.path$l1.total[keep.fits[[1]]], rule=2)$y
		cv.err.gamma[i,] <- approx(x=log10(path$gammas)[keep.fits[[(i+1)]]], y=cv.err[keep.fits[[(i+1)]]], xout=log10(orig.path$gammas)[keep.fits[[1]]], rule=2)$y

	}

	err.gamma <- rep(NA, n.gamma)
	err.gamma[keep.fits[[1]]] <- colSums(cv.err.gamma)
	err.se.gamma <- rep(NA, n.gamma)
	err.se.gamma[keep.fits[[1]]] <- sqrt(apply(cv.err.gamma, MARGIN=2, FUN=sd)/n.folds)
	cv.min.gamma <-  which.min(err.gamma)
	cv.1se.gamma.w <-  orig.path$gammas[ which(err.gamma < (err.gamma[cv.min.gamma] + err.se.gamma[cv.min.gamma]))]
	cv.1se.gamma <-  which(orig.path$gammas == max(cv.1se.gamma.w))

	err.l1 <- rep(NA, n.gamma)
	err.l1[keep.fits[[1]]]<- colSums(cv.err.l1)
	err.se.l1 <- rep(NA, n.gamma)
	err.se.l1[keep.fits[[1]]] <- sqrt(apply(cv.err.l1, MARGIN=2, FUN=sd)/n.folds)
	cv.min.l1 <-  which.min(err.l1)
	cv.1se.l1.w <-  orig.path$l1.total[which(err.l1 < (err.l1[cv.min.l1] + err.se.l1[cv.min.l1]))]
	cv.1se.l1 <- which(orig.path$l1.total==min(cv.1se.l1.w))

	return(list("sep.total"=orig.path$sep.total, "l1.total"=orig.path$l1.total,
			"cv.err.l1"=cv.err.l1, "cv.err.gamma"=cv.err.gamma, "err.gamma"=err.gamma, "err.se.gamma"=err.se.gamma, "err.l1"=err.l1, "err.se.l1"=err.se.l1,
			"cv.min.gamma"=cv.min.gamma, "cv.1se.gamma"=cv.1se.gamma,
			"cv.min.l1"=cv.min.l1, "cv.1se.l1"=cv.1se.l1,
			"n.test"=n.test, "gamma"=orig.path$gammas,
			"converged"=converged, "keep.fits"=keep.fits,
			"loss.type"=loss.type, "sep.list"=sep.list, "log.gamma.list"=log.gamma.list, "l1.list"=l1.list))

}

merge_close_regions <- function(z, margin=2){
	stopifnot(class(z)=="RangedData")
	starts <- c()
	ends <- c()
	scores <- c()
	ss <- unique(score(z))
	for(s in ss){
		q <- z[score(z) == s,]
		q <- reduce(q, min.gapwidth=margin)
		ends <- c(ends, end(q))
		starts <- c(starts, start(q))
		scores <- c(scores, rep(s, nrow(q)))
	}
	j <- order(starts)
	ir <- IRanges(starts[j], ends[j])
	rd <- RangedData(ir)
	rd$score <- scores[j]
	return(rd)
}

#Return table which can be written as bed file
#library(zoo)
get_separated_regions <- function(fit, which.window, chr="chr22", new.tol=NULL){
	require(IRanges)
	K <- dim(fit$fits)[2]
	stopifnot(K==3)

	if(!is.null(new.tol)){
		sep <- get_sep(fit$beta, new.tol)
		sep <- matrix(unlist(sep), ncol=3)
	}else{
		sep <- matrix(unlist(fit$sep), ncol=3)
	}
	#1-2, 1-3, 2-3
	if(all(sep==0)) return(0)

	S <- paste(sep[,1], sep[,2], sep[,3])
	z <- as(Rle(S), "RangedData")
	if(all( z$score == "0 0 0")) return(0)

	z <- z[ !score(z) == "0 0 0", ]
	#Merge over partition types
	z_merge <- reduce(z, min.gapwidth=2)
	z_merge <- z_merge[ width(z_merge) > 1, ]
	#Keep partition types separated
	z_sep <- merge_close_regions(z)
	z_sep <- z_sep[ width(z_sep) > 1, ]

	if(nrow(z_merge) ==0 & nrow(z_sep) ==0) return(0)

	res_tab_sep <- data.frame(matrix(nrow=nrow(z_sep), ncol=10))
	names(res_tab_sep) <- c("Chrom", "Start", "Stop", "Width", "Ncpg", "Window", "Partition", "AvgGap12", "AvgGap13", "AvgGap23")

	res_tab_merge <- data.frame(matrix(nrow=nrow(z_merge), ncol=9))
	names(res_tab_merge) <- c("Chrom", "Start", "Stop", "Width", "Ncpg", "Window", "AvgGap12", "AvgGap13", "AvgGap23")

	if(nrow(z_sep) > 0){
		res_tab_sep$Chrom <- chr
		res_tab_sep$Window <- which.window

		res_tab_sep$Start <- fit$pos[start(z_sep)]
		res_tab_sep$Stop <- fit$pos[end(z_sep)]
		res_tab_sep$Width <- res_tab_sep$Stop - res_tab_sep$Start
		res_tab_sep$Ncpg <- width(z_sep$ranges)

		for(i in 1:nrow(z_sep)){
			f <- fit$fits[start(z_sep)[i]:end(z_sep)[i],]
			p <- fit$pos[start(z_sep)[i]:(end(z_sep)[i])]

			f <- pmin(f, 1)
			f <- pmax(f, 0)

			k <- dim(f)[1]
			scoresum <- sum(as.numeric(unlist(strsplit(z_sep$score[i], split=" "))))

			auc <- apply(f, MARGIN=2, FUN=function(y){
				sum(diff(p)*rollmean(y,2))
			})
			if(z_sep$score[i] == "1 1 1" | scoresum== 1){

				g1 <- f[,1]-f[,2]
				res_tab_sep[i, "AvgGap12"]<-  sum(diff(p)*rollmean(abs(g1),2))/(max(p)-min(p))

				g2 <- f[,1]-f[,3]
				res_tab_sep[i, "AvgGap13"]<-  sum(diff(p)*rollmean(abs(g2),2))/(max(p)-min(p))

				g3 <- f[,2]-f[,3]
				res_tab_sep[i, "AvgGap23"]<-  sum(diff(p)*rollmean(abs(g3),2))/(max(p)-min(p))

				j <- order(auc)
				partition <- paste(c("(1)", "(2)", "(3)")[j], sep="", collapse="")
				if(scoresum==1) partition <- "Und"
				res_tab_sep$Partition[i] <- partition

			}else if(z_sep$score[i] == "1 1 0"){
				res_tab_sep$Partition[i] <- "(1)(2 3)"

				g <- f[,1]-f[,2]
				res_tab_sep[i, "AvgGap12"]<- res_tab_sep[i, "AvgGap13"]  <- sum(diff(p)*rollmean(abs(g),2))/(max(p)-min(p))

				if(auc[2] < auc[1]) res_tab_sep$Partition[i] <- "(2 3)(1)"

			}else if(z_sep$score[i] == "0 1 1"){
				res_tab_sep$Partition[i] <- "(1 2)(3)"

				g <- f[,1]-f[,3]
				res_tab_sep[i, "AvgGap13"]<-  res_tab_sep[i, "AvgGap23"] <- sum(diff(p)*rollmean(abs(g),2))/(max(p)-min(p))
				if(auc[3] < auc[1]) res_tab_sep$Partition[i] <- "(3)(1 2)"

			}else if(z_sep$score[i] == "1 0 1"){
				res_tab_sep$Partition[i] <- "(1 3)(2)"
				g <- f[,2]-f[,3]
				res_tab_sep[i, "AvgGap12"]<-  res_tab_sep[i, "AvgGap23"] <- sum(diff(p)*rollmean(abs(g),2))/(max(p)-min(p))

				if(auc[2] < auc[1]) res_tab_sep$Partition[i] <- "(2)(1 3)"

			}
		}
	}

	if(nrow(z_merge) > 0){
		res_tab_merge$Chrom <- chr
		res_tab_merge$Window <- which.window
		res_tab_merge$Start <- fit$pos[start(z_merge)]
		res_tab_merge$Stop <- fit$pos[end(z_merge)]
		res_tab_merge$Width <- res_tab_merge$Stop - res_tab_merge$Start
		res_tab_merge$Ncpg <- width(z_merge$ranges)
		for(i in 1:nrow(z_merge)){
			f <- fit$fits[start(z_merge)[i]:end(z_merge)[i],]
			p <- fit$pos[start(z_merge)[i]:end(z_merge)[i]]

			f <- pmin(f, 1)
			f <- pmax(f, 0)
			k <- dim(f)[1]

			g1 <- f[,1]-f[,2]
			res_tab_merge[i, "AvgGap12"]<-  sum(diff(p)*rollmean(abs(g1),2))/(max(p)-min(p))

			g2 <- f[,1]-f[,3]
			res_tab_merge[i, "AvgGap13"]<-  sum(diff(p)*rollmean(abs(g2),2))/(max(p)-min(p))

			g3 <- f[,2]-f[,3]
			res_tab_merge[i, "AvgGap23"]<-  sum(diff(p)*rollmean(abs(g3),2))/(max(p)-min(p))

		}
	}
	return(list("merged"=res_tab_merge, "separated"=res_tab_sep))

}


#library(plotrix)
plot_cv <- function(cv.obj, orig.path.file, meth.data, cv.file.name, fit.file.name, which.window=1, new.tol=NULL){
	k <- length(cv.obj$sep.list)
	png(cv.file.name, height=400*k, width=800)
	par(mfrow=c(k+1, 2))

	#Gamma v sep
	sep0 <- max(cv.obj$sep.list[[1]])
	x.ranges <- list()
	seprange <- c(0, sep0)
	#l1range <- range(unlist(lapply(1:k, FUN=function(x){ cv.obj$l1.list[[x]][cv.obj$keep.fits[[x]]]})))
	for(i in 1:k){
		#Try to find a reasonable x range
 		q <- lapply(1:length(cv.obj$log.gamma.list[[i]]), FUN=function(x){
				all(cv.obj$sep.list[[i]][cv.obj$log.gamma.list[[i]] > cv.obj$log.gamma.list[[i]][x]] < sep0)
				})
		q <- unlist(q)
		x.ranges[[i]] <- c(min(cv.obj$log.gamma.list[[i]][q]), max(cv.obj$log.gamma.list[[i]]))
		x.ranges[[i]][1] <- x.ranges[[i]][1]  - 0.2*(x.ranges[[i]][2]-x.ranges[[i]][1] )
		x.ranges[[i]][1] <- max(x.ranges[[i]][1], min(cv.obj$log.gamma.list[[i]][-1]))
		#cat(x.ranges[[i]], "\n")
		cols <- as.numeric(!cv.obj$keep.fits[[i]])+1
		plot(cv.obj$log.gamma.list[[i]], cv.obj$sep.list[[i]], xlim=x.ranges[[i]], ylim=seprange, ylab="Total Separated", xlab="Log Gamma", main=as.character(i-1), col=cols)

		l1range <- range(cv.obj$l1.list[[i]][ cv.obj$keep.fits[[i]] ])
		plot(cv.obj$log.gamma.list[[i]], cv.obj$l1.list[[i]], xlim=x.ranges[[i]],  ylim=l1range, ylab="L1 Distance", xlab="Log Gamma", main=as.character(i-1), col=cols)
	}

	#Cv err
	cat("CV err\n")
	l1range <- range(cv.obj$l1.list[[1]][ cv.obj$keep.fits[[1]] ])
	plotCI(x=cv.obj$l1.total, y=cv.obj$err.l1, uiw=cv.obj$err.se.l1, xlab="L1 Distance", ylab="CV Err", xlim=l1range)
	abline(v=cv.obj$l1.total[cv.obj$cv.min.l1], col="red", lty=2)
	abline(v=cv.obj$l1.total[cv.obj$cv.1se.l1], col="red", lty=1)
	dev.off()

	#cat("fit file: ", fit.file.name, "\n")
	png(fit.file.name, height=800, width=800)
	cat("Fit\n")
	path <- getobj(orig.path.file)
	fit <- path$JADE_fits[[cv.obj$cv.1se.l1]]

	R <- get_separated_regions(fit, which.window, new.tol=new.tol)
	if(length(R)==1){
		res.tab=NULL
	}else{
		 res.tab <- R$merged
	}

	dataplot <- plot_data(meth.data$y, meth.data$reads, pos=meth.data$pos)
	fits <- pmin(fit$fits, 1)
	fits <- pmax(fits, 0)
	fitplot <- plot_multi(fits=fits, wsize=100, pos=meth.data$pos, reads=meth.data$reads, ylim=c(0, 1), sep.tab=res.tab)
	par(mfrow=c(2, 1))
	multiplot(plotlist=list(dataplot, fitplot))
	dev.off()

	return(R)

}

