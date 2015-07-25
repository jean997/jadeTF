#'Estimate variance of a trendfiltering fit by bootstrapping
#'@param Y A matrix of complete (not average) data for one group.
#'@param READS A matrix of reads corresponding to \code{Y} (Optional). If \code{READS} is provided,
#'the elements in \code{Y} are assumed to be binomial observations and \code{sds} will be
#'calculated unless \code{calc.sds=FALSE}
#'@param sds A matrix of estimated standard deviations of \code{Y}. If missing
#'observations are assumed to have equal variance.
#'@param pos Vector of positions corresponding to the rows of \code{Y}
#'@param scale.pos
#'@param ord Order of trendfiltering. Can be 0, 1 or 2.
#'@param n.rep Number of bootstrap samples to take.
#'@return A list with three elements
#' \describe{
#'   \item{\code{all.fits}}{All of the fits from bootstrapped samples. p by n.rep matrx}
#'   \item{\code{avg.fit}}{Equivalent to \code{rowMeans(all.fits)}}
#'   \item{\code{var}}{Estimated variance.}
#' }
#'@export
bs_var_tf <- function(Y, lambda, sample.size,
										READS=NULL,  sds=NULL, calc.sds=TRUE,
										pos=NULL, scale.pos=NULL, ord=2,
										n.rep=100){

		p <- dim(Y)[1] #Number of sites
		k <- dim(Y)[2] #Number of samples

		#Scale positions
    if(!is.null(pos)){
      stopifnot(length(pos)==p)
    }else{
        pos <- 1:p
    }
		pos.give=pos
    if(!is.null(scale.pos)){
      R <-range(pos)
      pos <- scale.pos* ((pos-R[1])/(R[2]-R[1]))
    }

		if(is.null(sds)) sds <- rep(1, p)
		all.fits <- matrix(0, p, n.rep)
		i <- 1
		while(i <= n.rep){
			cat(i, "..")
			S <- sample(1:k, size=k, replace=TRUE)
			new.Y <- Y[, S]
			if(!is.null(READS)){
				new.READS <- READS[, S]
				y= rowSums(new.Y)/rowSums(new.READS)
				if(calc.sds){
				  yhat <- (rowSums(new.Y)+0.5)/(rowSums(new.READS)+1)
				  new.sds <- sqrt(yhat*(1-yhat)/rowSums(new.READS))
				}else{
				  new.sds <- sds
				}
			}else{
				y <- rowMeans(new.Y)
				new.sds <- sds
			}
			f <- fit_one(y, lambda, pos, new.sds, sample.size, ord=ord)
			all.fits[,i] <- f$fit
			i <- i+1
		}
		cat("\n")
		avg.fit <- rowMeans(all.fits)
		var.fit <- (1/(n.rep-1))*rowSums((all.fits-avg.fit)^2)
		return(list("all.fits"=all.fits, "avg"=avg.fit, "var"=var.fit))
}


