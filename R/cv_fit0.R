#' Run initial cross validation fits for JADE
#' @description Use an existing JADE fit at \eqn{\gamma = 0} to start cross validation fits with the same parameters.
#' These can be used as start computation of a series of solutions for \eqn{\gamma} between 0 and \eqn{\gamma_{max}}{\gamma_max}
#' using either \code{\link{jade_path_plan}} or \code{\link{jade_path}}.
#' @param orig.fit JADE fit using all the data or the name of a file containig a JADE fit.
#' @param n.folds Number of cross validation folds.
#' @param which.fold Vector of folds to fit. If This can be used to run cv fits in parallel.
#' If \code{which.fold=1:n.folds} or is left missing cv fits will be run sequentially.
#' @param data.file Optional file containing data. This is required if you want to bootstrap the variance.
#' @param n.rep.bs.var Number of reps to use in bootstrapping variance of the fit. If NULL varaince will
#' not be estimated.
#' @param save.prefix If provided each fit will be saved to file named save.prefix.fold.RData
#' @param return.objects If \code{TRUE} a list of jade objects will be returned. Otherwise
#' nothing will be returned.
#' @param lambda If \code{NULL} \eqn{\lambda} will be chosen by cross validation for each fold. Otherwise
#' \code{lambda} should be a list of length \code{which.fold} giving \eqn{lambda} for each fold.
#' @return A single jade object or list of jade objects of the same type as
#' \code{orig.fit} (either \code{jade_admm} or \code{jade_gd}) if \code{return.objects=TRUE}
#' @export
cv_fit0 <- function(orig.fit, n.folds=5, which.fold=1:n.folds, data.file=NULL,
                    n.rep.bs.var=NULL, save.prefix=NULL, return.objects=TRUE, lambda=NULL){

  if(class(orig.fit) == "character"){
    orig.fit.file <- orig.fit
    orig.fit <- getobj(orig.fit.file)
  }

  if(!is.null(lambda)){
    stopifnot(class(lambda)=="list")
    stopifnot(length(lambda)==length(which.fold))
  }
	if(!is.null(data.file) & !is.null(n.rep.bs.var)){
		R <- getobj(data.file)
		get.fit.var <- TRUE
	}else{
		get.fit.var <- FALSE
		if(!is.null(data.file) | !is.null(n.rep.bs.var)) cat("Warning: To bootstrap fit variance you must provide both data.file and n.rep.bs.var\n
		                                                     Variance will not be estimated.\n")
	}

	p <- dim(orig.fit$y)[1]
	K <- dim(orig.fit$y)[2]

	out <- list()
	out.ct <- 1
	folds <- matrix(0, nrow=p, ncol=K)
 	for(j in 1:K){
		non.missing <- which(!is.na(orig.fit$y[,j]))
		pj <- length(non.missing)
 		folds[non.missing,j] <- ((rep(1:n.folds, ceiling(pj/n.folds))[1:pj] + j ) %% n.folds ) + 1
 	}

	out <- list()
	out.ct <- 1
	for(i in which.fold){
	  cat(i, "..")
	  new.y <- orig.fit$y
	  test.idx <- matrix(0, nrow=p, ncol=K); test.idx[folds==i] <- 1
	  train.idx <- matrix(0, nrow=p, ncol=K); train.idx[folds==i & folds > 0] <- 1
    new.y[ folds==i] <- NA
    if(!is.null(lambda)) lambda.i <- lambda[[out.ct]]
        else lambda.i <- NULL
    if(orig.fit$algorithm=="admm"){
      fit <- jade_admm(y=new.y, gamma=0, pos=orig.fit$pos, scale.pos=orig.fit$scale.pos,
                       lambda=lambda.i, sample.size=orig.fit$sample.size, ord=orig.fit$ord,
                       sds=orig.fit$sds, tol=orig.fit$tol)
    }else if(orig.fit$algorithm=="gd"){
      fit <- jade_gd(y=new.y, gamma=0, pos=orig.fit$pos, scale.pos=orig.fit$scale.pos,
                       lambda1=lambda.i, lambda2=orig.fit$lambda2,
                       sample.size=orig.fit$sample.size, ord=orig.fit$ord,
                       sds=orig.fit$sds, sep.tol=orig.fit$tol)
    }

	  #Get fit.var
	  if(get.fit.var){
		  strt <- 1; stp <- 0
		  fit.var <- matrxi(nrow=p, ncol=K)
		  for(j in 1:K){
			  stp <- stp + orig.fit$sample.size[j]
			  if(is.null(R$READS)){
				  my.var <- bs_var_tf(R$Y[, strt:stp], lambda=fit0$lambda[j], sample.size=orig.fit$sample.size[j], positions=orig.fit$pos, scale.pos=orig.fit$scale.pos, ord=orig.fit$ord, n.rep=n.rep.bs.var, sds=orig.fit$sds[,j])
			  }else{
				  my.var <- bs_var_tf(R$Y[, strt:stp], lambda=fit0$lambda[j], sample.size=orig.fit$sample.size[j], positions=orig.fit$pos, scale.pos=orig.fit$scale.pos, ord=orig.fit$ord, n.rep=n.rep.bs.var, READS=R$READS[, strt:stp], sds=orig.fit$sds[,j])
			  }
			  fit.var[,j] <- my.var$var
			  strt <- strt + orig.fit$sample.size[j]
		  }
		  fit$fit.var <- fit.var
	  }
    if(!is.null(save.prefix)){
	    save(fit, file=paste(save.prefix, i, "RData", sep="."))
    }
    out[[out.ct]] <- fit
    out.ct <- out.ct + 1
	}
	if(length(which.fold)==1) out = out[[1]]
	if(return.objects) return(out)
	  else return()
}

