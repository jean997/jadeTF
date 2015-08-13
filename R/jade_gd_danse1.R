#' Fit JADE using gradient descent
#'
#' @description Solves the same problem solved by \code{\link{jade_admm}}.
#' This function also allows an additional L_1 penalty on the fitted values.
#'
#' @param y Data matrix of size p x K. May contain NA values but may not contain rows which are all NA.
#' @param gamma Fusion penalty.
#' @param pos Position vector of length p. If missing will use 1:p.
#' @param scale.pos An integer indicating to internally scale positions to range between 0 and \code{scale.pos}.
#' @param lambda1 Smoothing penalty vecor of length K.
#' If not provided, lambda will be chosen by cross validation.
#' @param lambda2 Parameter for L_1 penalty on fitted values. Defaults to 0.
#' @param sample.size Vector of sample sizes of length K.
#' If missing sample sizes are assumed to be 1.
#' @param ord Order of polynomial to fit. May be 0, 1, or 2.
#' @param sds Matrix of estimated standard deviations of size p x K.
#' These are the inverse of the diagonal elements of \eqn{A){i}}{A_i}.
#' Only the relative sizes of \code{sds} is important.
#' @param fit.var Matrix of size p x K of estimated variance of trendfiltering fits.
#' This will be used to construct the pairwise weight matrices \eqn{W}{W}.
#' Currently this is only supported for \eqn{K=2}.
#' \code{fit.var} can be estimated by bootstrapping. See \code{\link{bs_var}}.
#' @param var.wts If \code{fit.var} is not provided, the diagonal
#' elements of \eqn{W} may be specified here.
#' Since pairwise weights are currently only allowed for \eqn{K=2},
#' \code{var.wts} must be a vector of length p.
#' @param subset.wts This option can be used to obtain a de-biased fit with the
#' \eqn{\gamma} penalty only applied to pairs of points previously determined to be fused.
#' It should be a list of lists of the same format as the output of \code{\link{get_sep}}.
#' Elements are vectors of length p of 0s and 1s with 0 indicating
#' that the pair of points should not be penalized.
#' @param theta0,duals0 Starting values for \eqn{\theta} and the dual variable.
#'If a solution has been found for a nearby value
#' of \eqn{\gamma} using these values can improve convergence time.
#' If not provided the solution at \eqn{\gamma = 0} is used.
#' @param verbose Be chattier.
#' @param tol Tolerance for declaring points separated.
#' Separation can be recalculated with a different value of \code{tol} using \code{\link{get_sep}}.
#' @param max.it Maximum number of iterations.
#'
#' @return A \code{jade_tf} object. This really just a list with values including
#' \describe{
#'   \item{\code{fits}}{A p x K matrix of solutions.}
#'   \item{\code{n}}{Number of iterations to convergence}
#'  \item{\code{duals}} Dual variable
#'   \item{\code{sep}}{List of lists giving separation. See \code{\link{get_sep}}}
#' }
#' As well as all of the original parameters.
#' @export
jade_gd <- function(y, gamma, pos = NULL, scale.pos=NULL,
                           lambda1=NULL, lambda2=NULL, sample.size=NULL, ord=0,
													sds=NULL, fit.var=NULL, var.wts=NULL, subset.wts=NULL,
                          theta0=NULL, duals0=NULL, verbose=FALSE, sep.tol=1e-3,
													max.it=1000, thresh=1e-8, stepsize=NULL,
													eps=0, cv.metric=c("mse", "abs", "pois"), truncate.metric=100, shift=NULL, debug=FALSE){

  metric <- match.arg(cv.metric)

	stopifnot(ord %in% c(0, 1, 2))
	if(!is.null(var.wts) & !is.null(fit.var)) stop("Please provide only one of var.wts or fit.var")
	if(class(y)=="numeric"){
		p <- length(y)
		y <- matrix(y, nrow=p)
	}
	p <- dim(y)[1] #Number of sites
	K <- dim(y)[2] #Number of groups

	#Sample size
	if(is.null(sample.size)){
		cat("Assuming equal sample sizes are all 1 as sample.size is not provided.\n")
		sample.size <- rep(1, K)
	}

	#Data Summary
	if(verbose){
		cat("Groups:", K, "\n")
		cat("Sample sizes: ", sample.size, "\n")
		cat("Markers: ", p, "\n")
	}

	#Defaults
	#Standard deviations - A matrix
	if(is.null(sds)) sds <- matrix(1, p, K)

	#Fit var
	if(!is.null(fit.var)){
	  stopifnot(dim(fit.var) == c(p, K))
	  var.wts = wts_from_var(fit.var)
	}else if(!is.null(var.wts)){
	  stopifnot(length(var.wts)==(K-1))
	  for(j in 1:(K-1)) stopfinot(length(var.wts[[j]])==(K-j))
	}else{
	  var.wts = default_wts(p, K)
	}

	stopifnot(length(sample.size)==K)
	#Subset Weights
	if(!is.null(subset.wts)){
		if(K > 2) cat("Warning: subset.wts only used for 2 groups\n")
		stopifnot(class(subset.wts)=="list")
	}else{
		subset.wts=default_wts(p, K)
	}

	#Scale positions
	if(!is.null(pos)){
		stopifnot(length(pos)==p)
		pos.given <- pos
	}else{
		pos <- 1:p
		pos.given <- pos
	}
	if(!is.null(scale.pos)){
		R <-range(pos)
		pos<- scale.pos* ((pos-R[1])/(R[2]-R[1]))
	}

  ERRS <- 0 #Count errors


	if(is.null(lambda2)) lambda2 <- rep(0, K)


  ###Choose initial fits, cv lambda if necessary

  if(verbose) cat("Fitting at max value of gamma.\n")
  if(verbose & is.null(lambda1)) cat("Lambda1 will by chosen by cross validation.\n")
  theta.max <- fit_gammamax(y,  lambda1, pos, sample.size, sds, ord,
                      lambda2=lambda2, metric=metric, truncate.metric=truncate.metric, shift=shift)
  if(is.null(lambda1)) lambda1 <- theta.max$lambda
  theta.max <- theta.max$fit
  theta.min <- fit_gamma0(y,  lambda1,  pos, sample.size, sds, ord,
                      lambda2=lambda2)
  theta.min <- theta.min$fit
  if(!is.null(theta0)){
    theta <- theta0
  }else{
    theta <- theta.min
  }


	if(K==1 | gamma ==0){
		RETURN <- list("fits"=theta.min, "fit.max"=theta.max,
			               "y"=y, "sample.size"=sample.size, "fit.var"=fit.var,
			               "sds"=sds, "pos"=pos.given, "scale.pos"=scale.pos,
			               "lambda1"=lambda1, "lambda2"=lambda2, "gamma"=gamma, "ord"=ord,
			               "thresh"=thresh, "tol"=sep.tol, "subset.wts"=subset.wts, algorithm="gd")

		return(RETURN)
	}

  #Set stepsize
	if(is.null(stepsize)){
			stepsize <- min(0.5*min(sample.size))
	}

	#upper and lower bounds
	UB <- c(); LB <- c()

	#Initalize dual variables
	#Duals get initialized on the boundary if initial values not given
	if(is.null(duals0)){
		duals <-  starting_duals(theta, gamma, var.wts, sample.size, subset.wts)
	}else{
		duals <- constrain_duals(duals0, gamma, var.wts, sample.size, subset.wts)
	}

	#Iterate
  converged <- FALSE
	done <- FALSE
	iter.ct <- 1
	while(!done){
		lb0 <-dual_fct(y, theta, duals, lambda1, lambda2, sample.size, pos,  ord, sds)
		if(verbose) cat("lb0", lb0, "\n")
		#Update Theta
		tprev <- theta
		theta <- theta_update_tf(y, duals, sds, lambda1, lambda2, sample.size, pos, ord)
		#Check obj and dual values
		upper.bound <- obj_fct(y, theta, lambda1, lambda2, gamma, sample.size,
		                       subset.wts, sds, var.wts, pos, ord)
		lb1 <-dual_fct(y, theta, duals, lambda1, lambda2, sample.size, pos, ord, sds)
		if(verbose) cat("lb1", lb1, "\n")

    ##error##
		if((lb1-lb0)>1e-8){
			cat("Error! Theta update didn't reduce dual function: ", lb1-lb0, "\n")
			ERRS <- ERRS+1
      if(verbose){cat("Error..", ERRS, "\n")}
			if(debug){
			  R = list("y"=y, "tprev"=tprev, "theta"=theta, "duals"=duals, "sds"=sds, "lambda1"=lambda1,
			           "lambda2"=lambda2, "sample.size"=sample.size, "ord"=ord, "pos"=pos)
			  return(R)
			}
		}
			######
		UB <- c(UB, upper.bound)
		if(verbose) cat("upper bound: ", upper.bound, "\n")

		#Updadte duals
		for(j in 1:(K-1)){
			for(i in (j+1):K){
				u <- duals[[j]][[i-j]]
				#Gradient step
				#plus bigger-smaller
				u <- u + stepsize*(theta[,i]-theta[,j] - eps*u)
				#Box constraint
				w <- var.wts[[j]][[i-j]]
				u <- sign(u)*pmin(abs(u), gamma*w)
				#Weights
				wt <- subset.wts[[j]][[i-j]]
				u[wt==0] <- 0
        duals[[j]][[i-j]] <-u
			}
		}

		lower.bound <-dual_fct(y, theta, duals, lambda1, lambda2, sample.size, pos, ord, sds)
		LB <- c(LB, lower.bound)
		if(verbose) cat("lower bound: ", lower.bound, "\n")

		test1 <- (upper.bound - lower.bound)/upper.bound
		test2 <- max(abs(theta-tprev))
		if(verbose) cat("tests: ", test1, test2, "\n")
		stopifnot(test1 >= -1e-10)

		iter.ct <- iter.ct+1

		#if(test < thresh | test.1 < thresh )
		if(test1 < thresh & test2 < thresh){
			converged <- TRUE
			done <- TRUE
		}else if(iter.ct >= max.it){
		  done <- TRUE
		}
		if(verbose) cat(iter.ct, " ")
	}
	if(verbose) cat("\n")

	#if(penalty=="TF") cat("\n" )
	RETURN <- list("fits"=theta, "n"=iter.ct, "duals"=duals, "errors"=ERRS,"UB"=UB, "LB"=LB,
	               "y"=y, "sample.size"=sample.size, "sds"=sds,
	               "subset.wts"=subset.wts, "fit.var"=fit.var, "var.wts"=var.wts,
	               "pos"=pos.given, "scale.pos"=scale.pos, "lambda1"=lambda1, "lambda2"=lambda2,
	               "gamma"=gamma, "ord"=ord, "converged"=converged, "thresh"=thresh, "tol"=sep.tol,
							  "eps"=eps)
	sep <- sep_gd(RETURN, sep.tol)
	RETURN$sep <- sep
	return(RETURN)
}


#Minimize dual wrt thetas
#For each j solve
###
# n/2 || Ay - A\theta ||_2_^2 +  u\theta + lambda_1 || D\theta ||_1 + lambda_2 || \theta ||_1
###
# 1/2 (\theta^T A^T A \theta - (2y^TA^TA\theta - 2u/n)\theta   )
# 1/2 || (A^Ty - A^-1 u/n) - A\theta ||_2_^2 + lambda_1/n || D\theta ||_1 + lambda_2/n || \theta ||_1
theta_update_tf <- function(y, duals, sds, lambda1, lambda2, sample.size, pos, ord){
	K <- dim(y)[2]
	p <- dim(y)[1]

	fits <- matrix(0, nrow=p, ncol=K)

	#cat(p, K, "\n")
	weights <- 1/sds
	#weights[is.na(y)] <- 0

	for(j in 1:K){
		u <- rep(0, p)
		if(j < K){
			for(k in (j+1):K){
			  #j is smaller so substract
				u <-  u - duals[[j]][[k-j]]
			}
		}
		if(j > 1){
			for(k in 1:(j-1)){
			  #j is larger so add
				u <-  u + duals[[k]][[j-k]]
			}
		}

		nm <- which(!is.na(y[,j]))
		new.y <- as.vector(weights[,j]*y[,j] -( u/(weights[,j]*sample.size[j]) ))
		#cat(length(new.y),  ord, "\n")
		new.lam1 <- lambda1[j]/sample.size[j]
		if(all(weights[nm,j]==1)){
			out <- genlasso:::trendfilter(y=new.y[nm], pos=pos[nm], ord=ord)
		}else{
		  out <-  trendfilter_weights(y=new.y[nm], pos=pos[nm], ord=ord, wts=weights[nm,j])
		}
		co.nm <- coef.genlasso(out, lambda=new.lam1, type="primal")$beta

		if(length(nm) < p){
				co = .Call("tf_predict_R",
                    sBeta = as.double(co.nm),
                    sX = as.double(pos[nm]),
                    sN = length(nm),
                    sK = as.integer(ord),
                    sX0 = as.double(pos),
                    sN0 = length(pos),
                    sNLambda = 1,
                    sFamily = 0,
                    sZeroTol = as.double(1e-11),
                    package = "jadeTF")
		}else{
		  co <- co.nm
		}

		if(lambda2[j] > 0){
		  co <- soft_threshold(co, lambda2[j]/(sample.size[j]))
		}

		fits[,j] <- co
	}
  return(fits)
}

