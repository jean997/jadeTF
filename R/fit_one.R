
#Fit one sample
#Missing points are imputed
#Minimize N/(2) || (y - \theta)/sds ||^2 + \lambda_1||D\theta||_1 + \lambda_2||\theta||_1
#Equivalent to 1/2 || w*(y - \theta) ||^2 + \lambda_1/N||D\theta||_1 + \lambda2/N||\theta||_1
fit_one <- function(y, lambda, pos, sds, sample.size, ord,
                    lambda2=0, metric=c("mse", "abs", "pois"), truncate.metric=Inf, shift=NULL){

  metric <- match.arg(metric)

  p <- length(y)
  nm <- which(!is.na(y))
  wts <- 1/sds[nm]
  equal.wts <- all(wts == wts[1])
  #If all the weights are equal solve
  #Minimize 1/2 || y - \theta ||^2 + (\lambda_1/(N*w^2)||D\theta||_1

  if(equal.wts){
    tfit.out <- genlasso::trendfilter(y=y[nm], pos=pos[nm], ord=ord)
  }else{
    tfit.out <- trendfilter_weights(y=y[nm], pos=pos[nm], wts=wts, ord=ord)
  }

  if(is.na(lambda)){
    cv <- cv_pred.genlasso(obj=tfit.out, n.folds = 5, mode = "predict",
                           lambda2=lambda2/sample.size, metric=metric, truncate.metric=truncate.metric,
                           shift=shift)
    l <- cv$lambda.1se #l = lambda_1/(N*w^2) or lambda_1/N
    lambda <- l*sample.size
    if(equal.wts) lambda <- lambda*(wts[1])^2
    cat(lambda, "\n")
  }else if(equal.wts){
    l <- lambda/(sample.size*wts[1]^2)
  }else{
    l <- lambda/sample.size
  }
  co <- coef.genlasso(tfit.out, lambda = l, type="primal")$beta

  if(lambda2 > 0){
    co <- soft_threshold(co, lambda2/(sample.size*(wts^2)))
  }


  if(any(is.na(y))){
    fit = .Call("tf_predict_R",
                sBeta = as.double(co),
                sX = as.double(tfit.out$pos),
                sN = length(tfit.out$y),
                sK = as.integer(tfit.out$ord),
                sX0 = as.double(pos),
                sN0 = length(pos),
                sNLambda = 1,
                sFamily = 0,
                sZeroTol = as.double(1e-11), package="jadeTF")
  }else{
    fit <- co
  }
  return(list("fit"=fit, "lambda"=lambda))
}

