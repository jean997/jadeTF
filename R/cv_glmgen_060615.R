
#Takes a genlasso object but uses projection code from glmgen package to impute masked sites
#We wanted to minimize 1/2 || w(y - \theta) ||^2 + \lambda_1||D\theta||_1 + \lambda_2 ||\theta||_1
#Before running genlasso we created
# new_y = y*w and then computed trendfilter(new_y, pos=pos, X=diag(w), ord=ord)
# for a particular value of lambda_1 and lambda_2, the solution to the original problem can be found
# by soft thresholding the appropriate collumn of obj$beta at (1/w^2)*lambda_2
# These values should be compared to new_y/w to test mse

cv.genlasso_glmgen <- function(obj, weights, k=5, zero_tol=1e-11, mode=c("glmgen", "approx"), lambda2=0, metric=c("mse", "abs", "pois"), shift=NULL, truncate.metric=Inf){
	mode <- match.arg(mode)
	metric <- match.arg(metric)
	stopifnot("genlasso" %in% class(obj))

	p <- length(obj$y)
	folds <- c(0, rep(1:k, ceiling(p/k))[1:(p-2)], 0)
	cvall <- matrix(0, ncol=length(obj$lambda), nrow=k)
	family_cd = 0

	true_y <- obj$y/weights

	if(is.null(obj$pos)) obj$pos = 1:p

	for(i in 1:k){
			cat(i, "..")
			otr <- which(folds !=i)
			ote <- which(folds ==i)
			ytr <- obj$y[otr]

			xtr <- obj$pos[otr]
			xte <- obj$pos[ote]
			if(all(weights==1)){
				cvout <-  genlasso:::trendfilter(y=ytr, pos=xtr, ord=obj$ord)
			}else{
				A <- diag(weights[otr])
				cvout <-  genlasso:::trendfilter(y=ytr, pos=xtr, ord=obj$ord, X=A)
			}
			co = coef.genlasso(cvout, obj$lambda, type="primal")$beta
		#Soft threshold
		if(lambda2 > 0){
			co <- apply(co, MARGIN=2, FUN=function(x, l, weights){
			soft_threshold(x, (1/weights^2)*l)}, l=lambda2, weights=weights[otr])
		}
		X <- rbind(obj$lambda, co)
		#Predict
		if(mode=="approx"){
			cvpred <- apply(X, MARGIN=2, FUN=function(x, object, x.new){
							c <- x[-1]
							l <- x[1]
							y.new <- approx(x=object$pos, y=c, xout=x.new)$y
							return(y.new)},  object =cvout, x.new=xte)
		}else if(mode == "glmgen"){
			cvpred <- apply(X, MARGIN=2, FUN=function(x, object, x.new){
							c <- x[-1]
							l <- x[1]
							z = .Call("tf_predict_R",
                    sBeta = as.double(c),
                    sX = as.double(object$pos),
                    sN = length(object$y),
                    sK = as.integer(object$ord),
                    sX0 = as.double(x.new),
                    sN0 = length(x.new),
                    sNLambda = length(l),
                    sFamily = family_cd,
                    sZeroTol = as.double(zero_tol),
                    package = "jadeTF")
							return(z)}, object =cvout, x.new=xte)
		}
		if(metric == "mse"){
			cv_loss <-(true_y[ote]-cvpred)^2
			#truncate.metric <- truncate.metric^2
		}else if(metric=="abs"){
			cv_loss <-abs(true_y[ote]-cvpred)
		}else if(metric=="pois"){
			if(is.null(shift) & any(true_y < 0)) shift <- -1*min(true_y)
				else if(is.null(shift)) shift <- 0
			cv_loss <--(true_y[ote]+shift)*log(cvpred+shift)+cvpred+shift
		}
		cv_loss <- pmin(cv_loss, truncate.metric)
		cvall[i,] <- colMeans(cv_loss)
	}
	cat("\n")
	cverr <- colMeans(cvall)
	cvse <- apply(cvall, 2, sd)/sqrt(k)

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

