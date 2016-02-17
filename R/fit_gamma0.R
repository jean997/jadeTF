
#Fit JADE at gamma=0
fit_gamma0 <- function(y, lambda, pos, sample.size, sds, ord,
                       lambda2=NULL, metric=c("mse", "abs", "pois")){
  p <- dim(y)[1]
  K <- dim(y)[2]

  stopifnot(length(sample.size)==K)
  if(is.null(sds)) sds <- matrix(1, p, K)

  if(is.null(lambda)){
    lambda <- rep(NA, K)
  }else{
    stopifnot(length(lambda)==K)
  }

  if(is.null(lambda2)){
    lambda2= rep(0, K)
  }else{
    stopifnot(length(lambda2)==K)
  }

  fit <- matrix(0, p, K)
  for(j in 1:K){
    f <- fit_one(y[,j], lambda[j], pos, sds[,j], sample.size[j], ord=ord,
                 lambda2=lambda2[j], metric=metric)
    fit[,j] <- f$fit
    if(is.na(lambda[j])){
      lambda[j] <- f$lambda
    }
  }
  return(list("fit"=fit, "lambda"=lambda))
}
