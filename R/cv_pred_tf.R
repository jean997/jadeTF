#Cross validation for a trendfilter object with possible weights - can be produced using
  #genlasso::trendfilter() with X=NULL
  #or produced with jadeTF::trendfilter_weights

cv_pred.genlasso <- function(obj, n.folds=5, mode=c("predict", "approx"),
                             lambda2=0, metric=c("mse", "abs", "pois"), truncate.metric=Inf,
                             shift=NULL, zero.tol=1e-11){
	mode <- match.arg(mode)
	metric <- match.arg(metric)
	stopifnot("genlasso" %in% class(obj))

	p <- length(obj$y)
	nL <- length(obj$lambda)
	folds <- c(0, rep(1:n.folds, ceiling(p/n.folds))[1:(p-2)], 0)
	avg.test.loss <- matrix(0, ncol=nL, nrow=n.folds)

	if(is.null(obj$pos)) obj$pos = 1:p
  if(is.null(obj$weights)) wts = rep(1, p)
	  else wts = obj$weights
	for(i in 1:n.folds){
    cat(i, "..")
		otr <- which(folds !=i)
		ote <- which(folds ==i)

		ytr <- obj$y[otr]
    xtr <- obj$pos[otr]
    wtr <- wts[otr]

    xte <- obj$pos[ote]
    wte <- wts[ote]

		if(all(wts==1)){
			out.train <-  genlasso:::trendfilter(y=ytr, pos=xtr, ord=obj$ord)
		}else{
			out.train <-  trendfilter_weights(y=ytr, pos=xtr, ord=obj$ord, wts=wtr)
		}

    #matrix co.train is length(otr) x nL
		co.train = coef.genlasso(out.train, obj$lambda, type="primal")$beta

		#Soft threshold
		if(lambda2 > 0){
			co.train <- apply(co.train, MARGIN=2, FUN=function(x, l, wts){
			              soft_threshold(x, (1/wts^2)*l)},
			            l=lambda2, wts=wtr)
		}

		#Predict
		#co.test is length(ote) x nL
		if(mode=="approx"){
			co.test <- apply(co.train, MARGIN=2, FUN=function(c){
							y.new <- approx(x=xtr, y=c, xout=xte)$y
							return(y.new)})
		}else if(mode == "predict"){
			co.test <- apply(co.train, MARGIN=2, FUN=function(c){
							z = .Call("tf_predict_R",
                    sBeta = as.double(c),
                    sX = as.double(xtr),
                    sN = length(xtr),
                    sK = as.integer(obj$ord),
                    sX0 = as.double(xte),
                    sN0 = length(xte),
                    sNLambda = 1,
                    sFamily = 0,
                    sZeroTol = as.double(zero.tol),
                    package = "jadeTF")
							return(z)})
		}
		if(metric == "mse"){
			test.loss <-((obj$y[ote]-co.test)*wte)^2
		}else if(metric=="abs"){
			test.loss <-abs(obj$y[ote]-co.test)
		}else if(metric=="pois"){
			if(is.null(shift) & any(obj$y < 0)) shift <- -1*min(obj$y)
				else if(is.null(shift)) shift <- 0
			test.loss <- -(obj$y[ote]+shift)*log(co.test+shift)+co.test+shift
		}
		test.loss <- pmin(test.loss, truncate.metric)
		#Average over test points
		avg.test.loss[i,] <- colMeans(test.loss)
	}
	cat("\n")

	#Average over folds
	cverr <- colMeans(avg.test.loss)
	cvse <- apply(avg.test.loss, 2, sd)/sqrt(n.folds)

  names(cverr) = names(cvse) = round(obj$lambda,3)
  i0 = which.min(cverr)
  lam.min = obj$lambda[i0]
  lam.1se = max(obj$lambda[cverr<=cverr[i0]+cvse[i0]])
  i.min = which(obj$lambda==lam.min)
  i.1se = which(obj$lambda==lam.1se)

  out = list(err=cverr,se=cvse,lambda=obj$lambda,
      lambda.min=lam.min,lambda.1se=lam.1se,i.min=i.min,i.1se=i.1se)

	return(out)
}

