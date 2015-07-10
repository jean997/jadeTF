#' @import  clusterpathRcpp
#' @import  Matrix


#y should be a p by K matrix
#y can have NA values
#sds if provided should be p by K. Estimated standard deviation of each observation
#fit.var can be obtained from bootstrapping. aslo p by K.
	#Bootstrap estimate of the variance of the single group smooth. Only used if K=2
jade_admm <- function(y, gamma, pos=NULL, lambda=NULL, sample.size=NULL,
													sds=NULL, ord=2, fit.var=NULL, var.wts=NULL,
                          scale.pos=NULL, theta0=NULL, subset.wts=NULL,
													verbose=FALSE, tol=0.001,
													rho.alpha=NULL, rho.beta=1, max.it=1000,
													tau.decr=2, tau.incr=2, mu=10, e.rel=1e-4, e.abs=1e-8, adjust.rho.alpha=FALSE){

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
	if(!is.null(fit.var) | !is.null(var.wts)){
		if(K > 2) cat("Warning: Pairwise weight matrix W is only used for 2 groups in this implementation.\n")
	}
	#W matrix
	if(K==2){
		if(is.null(fit.var)) var.wts <- rep(1, p)
		  else var.wts = sqrt(fit.var[,1] + fit.var[,2])
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

	###Choose initial fits, cv lambda if necessary

	if(verbose) cat("Fitting at max value of gamma.\n")
	if(verbose & is.null(lambda)) cat("Lambda will by chosen by cross validation.\n")
	theta.max <- fit_gammamax(y=y,  lambda=lambda, pos=pos, sample.size=sample.size, sds=sds, ord=ord)
	if(is.null(lambda)) lambda <- theta.max$lambda
	theta.max <- theta.max$fit

	if(verbose) cat("Fitting at gamma=0\n")
	theta.min <- fit_gamma0(y=y,  lambda=lambda, pos=pos, sample.size=sample.size, sds=sds, ord=ord)
	theta.min <- theta.min$fit
	if(!is.null(theta0)){
		theta <- theta0
	}else{
		theta <- theta.min
	}

	if(K==1 | gamma ==0){
	  RETURN <- list("fits"=theta.min, "fit,max"=theta.max,
	                 "y"=y, "sample.size"=sample.size, "fit.var"=fit.var,
	                 "sds"=sds, "pos"=pos.given, "scale.pos"=scale.pos,
	                 "lambda"=lambda, "gamma"=gamma, "ord"=ord,
	                 "tol"=tol, "subset.wts"=subset.wts)

	  return(RETURN)
	}

	#Default starting value for step size - must be after lambda is chosen
	if(is.null(rho.alpha)) rho.alpha <- lambda*(((max(pos)-min(pos))/p)^(ord-1))

  #Prepare for ADMM algorithm
  #Build \tilde{D}^{k, x}
  if(ord==0){
    D <- diag(rep(1, p))
  }else{
	  D <- getDtfPosSparse(p, ord=ord-1, pos=pos)
		pos.wts <- ord/(pos[(ord+1):p]-pos[1:(p-ord)])
		D <- diag(pos.wts)%*%D
	}
	DtD <- crossprod(D)
	alpha.size <- dim(D)[1]

	AtA.diag <- get_AtA_diag(y, sds)

	AtAy <-  get_AtAy(y, sds)

	#Intitialize
	alpha <- D%*% theta
	u.alpha <- -D%*%theta +alpha

	beta <- theta
	u.beta <-  beta-theta

	#Inverses for theta_update
	theta.upd.qr.list <- list()
	for(j in 1:K){
		inv <- as(rho.alpha[j]*DtD + diag(rep(rho.beta, p) + (sample.size[j]*AtA.diag[,j])), "dgCMatrix")
		theta.upd.qr.list[[j]] <- qr(inv)
	}

	iter <- 0
	converged <- FALSE
	done <- FALSE

	#Start
	while(!done){
		theta.old <- theta
		beta.old <- beta
		alpha.old <- alpha

		#Update all varialbes
		theta <- theta_update(theta.upd.qr.list, AtAy, sample.size, AtA.diag, D, DtD, rho.alpha, rho.beta, alpha, u.alpha, beta, u.beta)
		alpha <- alpha_update(alpha, D, theta, u.alpha, lambda, rho.alpha)
		if(K==2) beta <- beta_update_k2(theta, u.beta, gamma, rho.beta, subset.wts, var.wts)
		  else beta <- beta_update(theta, u.beta, gamma, rho.beta)
		u.alpha <- u_alpha_update(u.alpha, theta, D, alpha)
		u.beta <- u_beta_update(u.beta, theta, beta)
		iter <- iter+1


		#Adjust rho.beta
		#Only the beta parts of residuals
		primal.resid.norm.beta <- sqrt(sum((beta-theta)^2))
		dual.resid.norm.beta <- sqrt(sum((rho.beta*(beta-beta.old))^2))
		if(primal.resid.norm.beta/dual.resid.norm.beta > mu){
			rho.beta <- tau.incr*rho.beta
			u.beta <- u.beta/tau.incr
			for(j in 1:K){
				theta.upd.qr.list[[j]] <- 0
			}
			if(verbose) cat("Changing rho.beta ", rho.beta, "\n")
		}else if(dual.resid.norm.beta/primal.resid.norm.beta > mu){
			rho.beta <- rho.beta/tau.decr
			u.beta <- u.beta*tau.decr
			for(j in 1:K){
				theta.upd.qr.list[[j]] <- 0
			}
			if(verbose) cat("Changing rho.beta ", rho.beta, "\n")
		}
		#if( iter %% 100 == 1){
		#	if(verbose) cat(iter, " ", primal.resid.norm.beta, " ", dual.resid.norm.beta, "\n")
		#}

		#Adjust rho.alphas - there is a different stepsize for each group
		primal.resid.norm.alpha <- sqrt(colSums((D%*%theta-alpha)^2))
		dual.resid.norm.alpha <- sqrt(rowSums((rho.alpha*(t(alpha-alpha.old)%*%D))^2))
		if(adjust.rho.alpha){
		  if(any(primal.resid.norm.alpha/dual.resid.norm.alpha > mu)){
		    idx <- which(primal.resid.norm.alpha/dual.resid.norm.alpha > mu)
		    for(j in idx){
		      rho.alpha[j] <- tau.incr*rho.alpha[j]
          u.alpha[,j] <- u.alpha[,j]/tau.incr
          theta.upd.qr.list[[j]] <- 0
          if(verbose) cat("Changing rho.alpha ",j, " ", rho.alpha[j], "\n")
		    }
		  }
		  if(any(dual.resid.norm.alpha/primal.resid.norm.alpha > mu)){
		  	idx <- which(dual.resid.norm.alpha/primal.resid.norm.alpha > mu)
			  for(j in idx){
				  rho.alpha[j] <- rho.alpha[j]/tau.decr
          u.alpha[,j] <- u.alpha[,j]*tau.decr
          theta.upd.qr.list[[j]] <- 0
          if(verbose) cat("Changing rho.alpha ",j, " ", rho.alpha[j], "\n")
			  }
		  }
		}

		primal.resid.norm.alpha <- sum(primal.resid.norm.alpha)
		dual.resid.norm.alpha <- sum(dual.resid.norm.alpha)

		if( iter %% 100 == 1){
			if(verbose) cat(iter, " ", primal.resid.norm.alpha, " ", dual.resid.norm.alpha, "\n")
		}


		#Calculate stopping criteria
		dual.resid.norm <- dual.resid.norm.alpha  + dual.resid.norm.beta
		primal.resid.norm <- primal.resid.norm.alpha  + primal.resid.norm.beta

		rel.cri.pri <- max( sqrt(sum((D%*%theta)^2)), sqrt(sum(alpha^2)) + sqrt(sum(beta^2)))
		rel.cri.dual <- sqrt(sum((t(D)%*%u.alpha)^2)) + sqrt(sum(u.beta^2))
		e.dual <- sqrt(K*(alpha.size + p))*e.abs + e.rel*rel.cri.dual
		e.primal <- sqrt(K*(alpha.size + p))*e.abs  + e.rel*rel.cri.pri
		if( primal.resid.norm < e.primal & dual.resid.norm < e.dual){
			converged <- TRUE
			done <- TRUE
		}else if(iter > max.it){
			done <- TRUE
		}
		#if( iter %% 100 == 1){
		#	if(verbose) cat(iter, " ", primal.resid.norm, " ", dual.resid.norm, "\n", e.primal, e.dual, "\n")
		#}
	}
	sep <- get_sep(beta, tol)

	RETURN <- list("fits"=theta, "n"=iter,  "D"=D, "alpha"=alpha, "beta"=beta, "rho.alpha"=rho.alpha, "rho.beta"=rho.beta,
							"y"=y, "sample.size"=sample.size, "u.beta"=u.beta, "u.alpha"=u.alpha,
							"sds"=sds, "fit.var"=fit.var, "pos"=pos.given, "scale.pos"=scale.pos,
							"lambda"=lambda, "gamma"=gamma, "ord"=ord, "sep"=sep, "tol"=tol,
							"likelihood"="Normal", "subset.wts"=subset.wts, "converged"=converged)
								#"t_theta"=t_theta, "t_beta"=t_beta, "t_alpha"=t_alpha, "t_ua"=t_ua, "t_ub"=t_ub)

	return(RETURN)
}


#Upddate functions
theta_update <- function(theta.upd.qr.list, AtAy, sample.size, AtA.diag, D, DtD, rho.alpha, rho.beta, alpha, u.alpha, beta, u.beta){
	K <- dim(AtAy)[2]
	p <- dim(AtAy)[1]
	new.theta <- matrix(0, p,K)
	for(j in 1:K){
		if(class(theta.upd.qr.list[[j]])=="numeric"){
			inv <- as(rho.alpha[j]*DtD + diag(rep(rho.beta, p) + (sample.size[j]*AtA.diag[,j])), "dgCMatrix")
			theta.upd.qr.list[[j]] <- qr(inv)
		}
		new.resp <- sample.size[j]*AtAy[,j] + rho.alpha[j]*t(D)%*%(alpha[,j]-u.alpha[,j]) + rho.beta*(beta[,j]-u.beta[,j])
		new.theta[,j] <- as.vector(qr.coef(theta.upd.qr.list[[j]],  new.resp))
	}
	return(new.theta)
}

alpha_update <- function(alpha, D, theta, u.alpha, lambda, rho.alpha){
	K <- dim(alpha)[2]
	alpha.size <- dim(alpha)[1]
	new.alpha <- matrix(0, alpha.size,K)
	sol <- rep(0, alpha.size)
	for(j in 1:K){
		new.resp <- as.vector(D%*%theta[,j] + u.alpha[,j])
		new.lam <- lambda[j]/rho.alpha[j]
		f <- .C("tf_dp_R", n = as.integer(alpha.size), y = as.double(new.resp), lam1 = as.double(new.lam), beta = as.double(sol))
    new.alpha[,j] <- f$beta
	}
	return(new.alpha)
}

beta_update <- function(theta, u.beta, gamma, rho.beta){
	X <- theta + u.beta #p by K
	cpath <- clusterpath.l1.id(t(X), LAPPLY=lapply)
	new.beta <- find_solution_clusterpath_cpp(cpath, lambda=gamma/rho.beta)
	return(new.beta)
}

beta_update_k2 <- function(theta, u.beta, gamma, rho.beta, subset.wts, var.wts){
	stopifnot(dim(theta)[2] ==2)
	wts <- subset.wts[[1]][[1]]
	var.wts <- var.wts[ wts==1]
	X <- theta + u.beta
	new.beta <- matrix(0, nrow=nrow(theta), ncol=ncol(theta))
	new.beta[wts==0,] <- X[wts==0,]
	X <- X[wts==1,]
	M <- (X[,1] + X[,2])/2
	new.y <-(X[,1]-X[,2])/2
	new.lam <- gamma/rho.beta
	D <- sign(new.y)*pmax(abs(new.y)-(var.wts*new.lam), 0)
	new.beta[wts==1,] <- cbind(M +D , M-D)
	return(new.beta)
}

u_alpha_update <- function(u.alpha, theta, D, alpha){
	K <- dim(alpha)[2]
	alpha.size <- dim(alpha)[1]
	new.u.alpha <- u.alpha + (D%*%theta - alpha)
}

u_beta_update <- function(u.beta, theta, beta){
	K <- dim(beta)[2]
	p <- dim(beta)[1]
	new.u.beta <- u.beta + (theta - beta)
}
