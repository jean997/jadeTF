#' Fit JADE
#'
#' @description This is the main function in the jadeTF package. It fits JADE by minimizing
#' \deqn{\sum_{i=1}^{K} \Big( \frac{n_{i}}{2}\Vert A_{i}(y_{i} - \theta_{i}) \Vert_{2}^{2}  +
#' \lambda_{i}\Vert D^{k+1, x}\theta_{i}\Vert_{1}\Big) +
#' \gamma \sum_{i < l}\Vert W_{il} (\theta_{i} - \theta_{l}) \Vert_{1}}{
#' \sum_i=1^K ( n_i/2 || A_i y_i - \theta_i ||_2^2  +
#' \lambda_i || D \theta_i||_1) +
#' \gamma \sum_{j < l}|| W_jl (\theta_j - \theta_l) ||_1}
#'
#' @param y Data matrix of size p x K. May contain NA values but may not contain rows which are all NA.
#' @param gamma Fusion penalty.
#' @param pos Position vector of length p. If missing will use 1:p.
#' @param scale.pos An integer indicating to internally scale positions to range between 0 and \code{scale.pos}.
#' @param lambda Smoothing penalty vecor of length K.
#' If not provided, lambda will be chosen by cross validation.
#' @param sample.size Vector of sample sizes of length K.
#' If missing sample sizes are assumed to be 1.
#' @param ord Order of polynomial to fit. May be 0, 1, or 2.
#' @param sds Matrix of estimated standard deviations of size p x K.
#' These are the inverse of the diagonal elements of \eqn{A){i}}{A_i}.
#' Only the relative sizes of \code{sds} is important.
#' @param fit.var Matrix of size p x K of estimated variance of trendfiltering fits.
#' This will be used to construct the pairwise weight matrices \eqn{W}{W}.
#' Currently this is only supported for \eqn{K=2}.
#' \code{fit.var} can be estimated by bootstrapping. See \code{\link{bs_var_tf}}.
#' @param var.wts If \code{fit.var} is not provided, the diagonal
#' elements of \eqn{W} may be specified here.
#' Since pairwise weights are currently only allowed for \eqn{K=2},
#' \code{var.wts} must be a vector of length p.
#' @param subset.wts This option can be used to obtain a de-biased fit with the
#' \eqn{\gamma} penalty only applied to pairs of points previously determined to be fused.
#' It should be a list of lists of the same format as the output of \code{\link{get_sep}}.
#' Elements are vectors of length p of 0s and 1s with 0 indicating
#' that the pair of points should not be penalized.
#' @param theta0,u.alpha0,u.beta0 Starting values for \eqn{\theta} and the dual variables.
#'If a solution has been found for a nearby value
#' of \eqn{\gamma} using these values can improve convergence time.
#' If not provided the solution at \eqn{\gamma = 0} is used.
#' @param verbose Be chattier.
#' @param tol Tolerance for declaring points separated.
#' Separation can be recalculated with a different value of \code{tol} using \code{\link{get_sep}}.
#' @param max.it Maximum number of iterations.
#' @param rho.alpha,rho.beta ADMM step sizes.
#' \code{rho.alpha} has length K and \code{rho.beta} is a constant.
#' Change with caution.
#' @param tau.incr,tau.decr,mu Parameters for adjusting step size. Change with caution.
#' @param e.rel,e.abs Parameters for determining convergence. Change with caution.
#' @param adjust.rho.alpha Adaptively change rho.alpha.
#' This does not seem to help so the default is FALSE.
#'
#' @return A \code{jade_tf} object. This really just a list with values including
#' \describe{
#'   \item{\code{fits}}{A p x K matrix of solutions.}
#'   \item{\code{n}}{Number of iterations to convergence}
#'   \item{\code{beta},\code{alpha}}{See paper.}
#'   \item{\code{u.alpha}{u.beta}}{Dual variables. See paper.}
#'   \item{\code{sep}}{List of lists giving separation. See \code{\link{get_sep}}}
#' }
#' As well as all of the original parameters.
#' @export
jade_admm <- function(y, gamma, pos=NULL, scale.pos=NULL, lambda=NULL, sample.size=NULL, ord=2,
                      sds=NULL, fit.var=NULL, var.wts=NULL, subset.wts=NULL,
                      theta0=NULL, u.alpha0 = NULL, u.beta0=NULL,
                      verbose=FALSE, tol=0.001, max.it=1000,
                      rho.alpha=NULL, rho.beta=1, adjust.rho.alpha=TRUE,
                      tau.decr=2, tau.incr=2, mu=10, e.rel=1e-4, e.abs=1e-8){

  stopifnot(ord %in% c(0, 1, 2))
  if(!is.null(var.wts) & !is.null(fit.var)) stop("Please provide only one of var.wts or fit.var")
  if(class(y)=="numeric"){
    p <- length(y)
    y <- matrix(y, nrow=p)
    if(!is.null(sds)) sds <- matrix(sds, nrow=p)
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
  stopifnot(all(sds > 0, na.rm=TRUE))
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
    RETURN <- list("fits"=theta.min, "fit.max"=theta.max,
                   "y"=y, "sample.size"=sample.size, "fit.var"=fit.var,
                   "sds"=sds, "pos"=pos.given, "scale.pos"=scale.pos,
                   "lambda"=lambda, "gamma"=gamma, "ord"=ord,
                   "tol"=tol, "subset.wts"=subset.wts, algorithm="admm")

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
  if(is.null(u.alpha0)) u.alpha <- -D%*%theta +alpha
    else u.alpha <- u.alpha0

  beta <- theta
  if(is.null(u.beta0)) u.beta <-  matrix(0, nrow=p, ncol=K)
    else u.beta <- u.beta0

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

    #Adjust rho.alphas - there is a different stepsize for each group
    primal.resid.norm.alpha <- sqrt(colSums((D%*%theta-alpha)^2)) #Lenth K
    dual.resid.norm.alpha <- sqrt(rowSums((rho.alpha*(t(alpha-alpha.old)%*%D))^2))
    if(adjust.rho.alpha){
      #Increase rho alphas?
      if(any(primal.resid.norm.alpha> dual.resid.norm.alpha*mu)){
        idx <- which(primal.resid.norm.alpha/dual.resid.norm.alpha > mu)
        for(j in idx){
          rho.alpha[j] <- tau.incr*rho.alpha[j]
          u.alpha[,j] <- u.alpha[,j]/tau.incr
          theta.upd.qr.list[[j]] <- 0
          if(verbose) cat("Changing rho.alpha ",j, " ", rho.alpha[j], "\n")
        }
      }else if(any(dual.resid.norm.alpha > primal.resid.norm.alpha*mu)){
        #Decrease rho alphas?
        idx <- which(dual.resid.norm.alpha/primal.resid.norm.alpha > mu)
        for(j in idx){
          rho.alpha[j] <- rho.alpha[j]/tau.decr
          u.alpha[,j] <- u.alpha[,j]*tau.decr
          theta.upd.qr.list[[j]] <- 0
          if(verbose) cat("Changing rho.alpha ",j, " ", rho.alpha[j], "\n")
        }
      }
    }


    #Calculate stopping criteria
    dual.resid.norm <- sqrt( sum((rho.alpha*(t(alpha-alpha.old)%*%D))^2)  +  sum((rho.beta*(beta-beta.old))^2))
    primal.resid.norm <- sqrt(  sum((D%*%theta-alpha)^2)   +  sum((beta-theta)^2) )

    rel.cri.pri <- max( sqrt(sum((D%*%theta)^2) + sum(theta^2)), sqrt(sum(alpha^2) + sum(beta^2)))
    rel.cri.dual <- sqrt(sum((t(D)%*%u.alpha)^2) + sum(u.beta^2))
    e.dual <- sqrt(K*(alpha.size + p))*e.abs + e.rel*rel.cri.dual
    e.primal <- sqrt(K*(alpha.size + p))*e.abs  + e.rel*rel.cri.pri

    if(verbose & iter%%100==1){
      cat(primal.resid.norm, " ", e.primal, "\n")
      cat(dual.resid.norm, " ", e.dual, "\n")
      cat((primal.resid.norm < e.primal & dual.resid.norm < e.dual), "\n")
      cat(max(abs(theta.old-theta)), "\n")
    }
    if( primal.resid.norm < e.primal & dual.resid.norm < e.dual & iter > 1){
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

  RETURN <- list("fits"=theta, "sep"=sep, "n"=iter,  "D"=D,
                 "alpha"=alpha, "beta"=beta, "u.beta"=u.beta, "u.alpha"=u.alpha,
                 "rho.alpha"=rho.alpha, "rho.beta"=rho.beta,
                 "y"=y, "sample.size"=sample.size,
                 "sds"=sds, "fit.var"=fit.var, "pos"=pos.given, "scale.pos"=scale.pos,
                 "lambda"=lambda, "gamma"=gamma, "ord"=ord, "tol"=tol,
                 "subset.wts"=subset.wts, "converged"=converged, algorithm="admm")
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
