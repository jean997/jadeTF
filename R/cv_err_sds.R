#' Use cross validation paths to choose a value of gamma
#'
#' This function requires having a \code{jade_path} object for the original data
#' and for data with missing folds. See \code{\link{cv_fit0}}.
#'
#' @param orig.path Either a path object or a file containing a \code{jade_path} object.
#' @param cv.path.list Either a list of \code{jade_path} objects or a vector of file names
#' containing \code{jade_path} objects for the cross validation data sets.
#' @param use.converged.only Only use fits which have converged.
#' @param control.l1 Only use fits with \code{l1.total <= l1.total0}.
#' @return A list with elements
#' #' \describe{
#'  \item{\code{l1.total}}{A vector of length N where N is the number of fits in \code{orig.path} giving the total
#'  L1 distance between all pairs of profiles. Equivalent to \code{orig.path$l1.total}.}
#'  \item{\code{cv.err.l1}}{An \code{n.folds} by N matrix where N is the number of fits in \code{orig.path}.
#'  Each row gives the average corss validation error for each value of \code{l1.total}.}
#'  \item{\code{err.l1}}{A vector of length N giving average cross validation error over all folds.
#'  Equivalent to \code{colSums(cv.err.l1)} with missing values for fits that were discarded according to
#'  specification of \code{use.converged.only} and \code{control.l1}.}
#'  \item{\code{err.se.l1}}{Estimated standard error of \code{err.l1}.}
#'  \item{\code{cv.min.l1,cv.1se.l1}}{Indices of the fits in \code{orig.path} corresponding to the minimum and 1
#'  standard error rule cross validation error.}
#' }
#' @export
cv_err_wts <- function(orig.path, cv.path.list=NULL,
                       use.converged.only=TRUE, control.l1=TRUE){

  if(class(orig.path)=="character"){
    orig.path.file <- orig.path
    orig.path <- getobj(orig.path.file)
  }
  if(class(cv.path.list) == "character") files <- TRUE
    else files <- FALSE

  n.folds=length(cv.path.list)

	#Path with full data
	orig.y <- orig.path$JADE_fits[[1]]$y
	orig.sds <- orig.path$JADE_fits[[1]]$sds
	orig.na <- which(is.na(orig.y))
	K <- dim(orig.y)[2]
	p <- dim(orig.y)[1]
	n.gamma <- length(orig.path$JADE_fits)


	#Track convergence
	converged <- list(orig.path$converged)

	keep.fits <- list(rep(TRUE, n.gamma))
	if(use.converged.only){
		keep.fits[[1]][!converged[[1]]]<- FALSE
	}
	if(control.l1){
		keep.fits[[1]][ orig.path$l1.total > orig.path$l1.total[1]] <- FALSE
	}

	#CV by matching l1 distance
	cv.err.l1 <-  matrix(0, nrow=n.folds, ncol=sum(keep.fits[[1]]))

	log.gamma.list <- list(log10(orig.path$gammas))
	l1.list <- list(orig.path$l1.total)
	sep.list <- list(orig.path$sep.total)

	n.test <- c()

	#Collect info from each fold
	for(i in 1:n.folds){
	  if(files) path <- getobj(cv.path.list[i])
	    else path <- cv.path.list[[i]]

		cv.y <- path$JADE_fits[[1]]$y


		#convergence
		converged[[(i+1)]] <- path$converged

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

		cv.err <- unlist( lapply(path$JADE_fits, FUN=function(x, orig.y, orig.sds, test.idx){
					sum(( (x$fits[test.idx]-orig.y[test.idx])/orig.sds[test.idx] )^2)
					}, orig.y=orig.y, orig.sds=orig.sds, test.idx=test.idx))/n.test[i]

		cv.err.l1[i,] <- approx(x=path$l1.total[keep.fits[[(i+1)]]],
		                        y=cv.err[keep.fits[[(i+1)]]],
		                        xout=orig.path$l1.total[keep.fits[[1]]],
		                        rule=2)$y
	}


	err.l1 <- rep(NA, n.gamma)
	err.l1[keep.fits[[1]]]<- colMeans(cv.err.l1)
	err.se.l1 <- rep(NA, n.gamma)

	err.se.l1[keep.fits[[1]]] <- apply(cv.err.l1, MARGIN=2, FUN=sd)/sqrt(n.folds)
	cv.min.l1 <-  which.min(err.l1)
	if(length(cv.min.l1) > 1){
	  g <- orig.path$gammas[cv.min.l1]
	  cv.min.l1 <- cv.min.l1[ which.min(g)]
	}
	cv.1se.l1.w <-  orig.path$l1.total[which(err.l1 < (err.l1[cv.min.l1] + err.se.l1[cv.min.l1]))]
	cv.1se.l1 <- which(orig.path$l1.total==min(cv.1se.l1.w))
	if(length(cv.1se.l1) > 1){
	  g <- orig.path$gammas[cv.1se.l1]
	  cv.1se.l1 <- cv.1se.l1[ which.min(g)]
	}
	return(list("sep.total"=orig.path$sep.total, "l1.total"=orig.path$l1.total,
			"cv.err.l1"=cv.err.l1,  "err.l1"=err.l1, "err.se.l1"=err.se.l1,
			"cv.min.l1"=cv.min.l1, "cv.1se.l1"=cv.1se.l1,
			"n.test"=n.test, "gamma"=orig.path$gammas,
			"converged"=converged, "keep.fits"=keep.fits,
			"sep.list"=sep.list, "l1.list"=l1.list, "log.gamma.list"=log.gamma.list))

}
