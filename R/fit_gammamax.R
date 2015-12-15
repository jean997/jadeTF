#Fit JADE at gamma max
fit_gammamax <- function(y, lambda, pos, sample.size, sds, ord,
                         lambda2=NULL, metric=c("mse", "abs", "pois"), truncate.metric=Inf, shift=NULL){
  p <- dim(y)[1]
  K <- dim(y)[2]

  stopifnot(length(sample.size)==K)
  if(is.null(sds)) sds <- matrix(1, p, K)


  miss <- is.na(y)

  y[miss] <- 0
  sds[miss] <- 0
  ss <- matrix(rep(sample.size, p), byrow=TRUE, nrow=p)
  z <- ss/(sds^2); z[miss] <- 0
  new.sigma <- sqrt( 1/ rowSums(z) )
  new.y <- (y*ss)/(sds^2); new.y[miss] <- 0
  new.y <- rowSums(new.y)*(new.sigma^2)

  if(is.null(lambda)){
    new.lam <- NA
  }else{
    stopifnot(length(lambda)==K)
    new.lam <- sum(lambda)
  }

  if(is.null(lambda2)){
    new.lam2 <- 0
  }else{
    stopifnot(length(lambda2)==K)
    new.lam2 <- sum(lambda2)
  }
  fit <- fit_one(new.y, new.lam, pos, new.sigma, 1, ord,
                 lambda2=new.lam2, metric=metric, truncate.metric=truncate.metric, shift=shift)
  if(is.null(lambda)){
    lambda <- rep(fit$lambda/K, K)
  }

  fit <- matrix(rep(fit$fit, K), byrow=FALSE, nrow=p)
  return(list("fit"=fit, "lambda"=lambda))
}
