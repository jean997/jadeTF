#' @useDynLib jadeTF
#' @importFrom Rcpp sourceCpp
#' @import  genlasso


#Contains small utility functions

h <- function(x, p, pos, ord){
	D_tf <- getDtfPosSparse(n=p, ord=ord, pos=pos)
	return(sum(abs ( D_tf %*% x)))
}

#JADE Objective function for trend filtering problem
obj_fct <-  function(y, tfits, lambda, gamma, sample.size, subset.wts, fit_var, pos, ord){

		p <- dim(y)[1]
		K <- dim(y)[2]

		obj_value <- 0
		for(j in 1:K){
			obj_value <- obj_value + (sample.size[j]/2)*sum((y[,j]-tfits[,j])^2) + lambda[j]*h(tfits[,j], p, pos, ord)
			if(j==K) next
			for (i in (j+1):K){
				pen <- gamma*sum(abs(tfits[,j]-tfits[,i])*subset.wts[[j]][[i-j]]*(sqrt(fit_var[,i] + fit_var[,j])))
				obj_value <- obj_value + pen
			}
		}
		return(obj_value)
}

#Fit one sample
#Missing points are imputed
#Minimize 1/(2*N) || y - \theta ||^2 + \lambda_1||D\theta||_1
#Equivalent to 1/2 || y - \theta ||^2 + \lambda_1*N||D\theta||_1
fit_one <- function(y, lambda, pos, sds, sample.size, ord){
    p <- length(y)
		nm <- which(!is.na(y))
		wts <- 1/sds[nm]
		equal.wts <- all(wts == wts[1])
		#If all the weights are equal solve
		#Minimize 1/2 || y - \theta ||^2 + (\lambda_1*N)/w^2||D\theta||_1

		if(equal.wts){
			tfit.out <- genlasso:::trendfilter(y=y[nm], pos=pos[nm], ord=ord)
		}else{
		  tfit.out <- trendfilter_weights(y=y[nm], pos=pos[nm], wts=wts, ord=ord)
		}

    if(is.na(lambda)){
      cv <- cv_pred.genlasso(obj=tfit.out, n.folds = 5, mode = "predict")
      l <- cv$lambda.1se
      lambda <- l/sample.size
      if(equal.wts) lambda <- lambda*(weights[1])^2
      cat(lambda, "\n")
		}else if(equal.wts){
			l <- lambda/(weights[1]^2)
		}else{
			l <- lambda*sample.size
		}
		co <- coef.genlasso(tfit.out, lambda = l, type="primal")$beta
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
		              sZeroTol = as.double(1e-6), package="jadeTF")
		}else{
		  fit <- co
		}
    return(list("fit"=fit, "lambda"=lambda))
}


#Fit JADE at gamma=0
fit_gamma0 <- function(y, lambda, pos, sample.size, sds, ord){
		p <- dim(y)[1]
		K <- dim(y)[2]

		stopifnot(length(sample.size)==K)
		if(is.null(sds)) sds <- matrix(1, p, K)

		if(is.null(lambda)){
				lambda <- rep(NA, K)
		}else{
				stopifnot(length(lambda)==K)
		}

		fit <- matrix(0, p, K)
		for(j in 1:K){
				f <- fit_one(y[,j], lambda[j], pos, sds[,j], sample.size[j], ord=ord)
				fit[,j] <- f$fit
				if(is.na(lambda[j])){
						lambda[j] <- f$lambda #Accomodating sample.size occurs in fit_one
				}
		}
		return(list("fit"=fit, "lambda"=lambda))
}

#Fit JADE at gamma max
fit_gammamax <- function(y, lambda, pos, sample.size, sds, ord){
		p <- dim(y)[1]
		K <- dim(y)[2]

		stopifnot(length(sample.size)==K)
		if(is.null(sds)) sds <- matrix(1, p, K)


		missing <- is.na(y)

		y[missing] <- 0
		sds[missing] <- 0
		ss <- matrix(rep(sample.size, p), byrow=TRUE, nrow=p)
		z <- ss/(sds^2); z[missing] <- 0
		new_sigma <- sqrt( 1/ rowSums(z) )
		new_y <- (y*ss)/(sds^2); new_y[missing] <- 0
		new_y <- rowSums(new_y)*(new_sigma^2)
		if(is.null(lambda)){
			new_lam <- NA
		}else{
				stopifnot(length(lambda)==K)
				new_lam <- sum(lambda)
		}

		fit <- fit_one(new_y, new_lam, pos, new_sigma, 1, ord=ord)
		if(is.null(lambda)){
			l <- fit$lambda/sum(sample.size)
			lambda <- sample.size * l
		}

		fit <- matrix(rep(fit$fit, K), byrow=FALSE, nrow=p)
		return(list("fit"=fit, "lambda"=lambda))
}

#Default weights
default_wts <- function(p, K){
	wts <- list()
  for(j in 1:(K-1)){
    wts[[j]] <- list()
    for(i in (j+1):K){
      wts[[j]][[i-j]] <- rep(1, p)
    }
	}
	return(wts)
}

pairwise_wts <- function(subset.wts, fit_var, sample.size){
	K <- length(subset.wts)+1
	pw <- list()
	for(j in 1:(K-1)){
		pw[[j]] <- list()
		for(i in (j+1):K){
			pw[[j]][[i-j]] <- subset.wts[[j]][[i-j]]*sqrt( fit_var[i] + fit_var[j])
		}
	}
	return(pw)
}

soft_threshold <- function(vec,lam){
  # Soft threshold function
  #
  # ARGUMENTS
  #	vec	: vector that soft tresholding is applied
  #	lam	: non-negative scalar soft tresholding parameter.
  #
  # VALUES
  #	res	: resulting vector of the same size as vec
  #
  if ( length(lam)>1 &  length(lam)!=length(vec)){
    cat(length(vec), " ", length(lam), "\n")
    cat('\n ERROR: THE SIZE OF THE SECOND ARGUMENT SHOULD BE 1 or length of vec.\n')
    return ( 0 )
  }

  idx.1<-which(vec < -lam)
  idx.2<-which(vec > lam)
  res<-rep(0,length(vec))

  if(length(lam)==1){
    if ( length(idx.1)>0 ) res[idx.1]<- vec[idx.1]+lam
    if ( length(idx.2)>0 ) res[idx.2]<- vec[idx.2]-lam
  }else{
    if ( length(idx.1)>0 ) res[idx.1]<- vec[idx.1]+lam[idx.1]
    if ( length(idx.2)>0 ) res[idx.2]<- vec[idx.2]-lam[idx.2]
  }
  return( res )
}


#Helper functions
get_AtA_diag <- function(y, sds){
	K <- dim(y)[2]
	p <- dim(y)[1]
	AtA_diag <- matrix(0, p, K)
	for(j in 1:K){
		w <- rep(1, p)
		w[is.na(y[,j])] <- 0
		s <- sds[!is.na(y[,j]), j]
		w[!is.na(y[,j])] <- w[!is.na(y[,j])]/(s^2)
		AtA_diag[,j] <- w
	}
	return(AtA_diag)
}

get_AtAy <- function(y, sds){
	p <- dim(y)[1]
	K <- dim(y)[2]
	AtAy <- matrix(0, p, K)
	for(j in 1:K){
		idx <- which(!is.na(y[,j]))
		n.idx <- which(is.na(y[,j]))
		w <- rep(1, p)
    w[n.idx] <- 0
    w[idx] <- w[idx]/(sds[idx,j]^2)
		newy <- y[,j]
		newy[n.idx] <- 0
		AtAy[,j] <- newy*w
	}
	return(AtAy)
}


#Coppied from GWAS Tools
getobj <- function (Rdata){
    objname <- load(Rdata)
    if (length(objname) > 1) {
        warning(paste("Multiple objects stored in file", Rdata,
            "\nReturning only the first object"))
    }
    return(get(objname))
}






get_sep <- function(fits, tol){
	K <- dim(fits)[2]
	p <- dim(fits)[1]
	sep <- list()
	for(j in 1:(K-1)){
			sep[[j]] <- list()
			for(i in (j+1):K){
				u <- abs(fits[,i]- fits[,j])
				sep[[j]][[i-j]] <- rep(0, p)
				sep[[j]][[i-j]][u > tol] <- 1
			}
	}
	return(sep)
}

fused_from_sep <- function(sep){
	K <- length(sep)+1
	fused <- list()
	for(j in 1:(K-1)){
      fused[[j]] <- list()
      for(i in (j+1):K){
        fused[[j]][[i-j]] <- 1-sep[[j]][[i-j]]
      }
  }
  return(fused)
}

obj_fct_parts <-  function(fit.obj){

		y <- fit.obj$y
		fits <- fit.obj$fits
		sds <- fit.obj$sds
		pos <- fit.obj$pos
		ord <- fit.obj$ord
		p <- dim(y)[1]
		K <- dim(y)[2]

		alph_p <- dim(fit.obj$D)[1]
		D1 <- getDtf(n=alph_p, ord=0)
		D <- D1 %*% fit.obj$D

		miss <- which(is.na(y))
		sds[miss] <- Inf
		ssd <- sum(((y-fits)/sds)^2, na.rm=TRUE)/2
		tf_pen <- rep(0, K)
		cl_pen <- 0
		for(j in 1:K){
			tf_pen[j] <- fit.obj$lambda[j]*sum(abs(D%*%fits[,j]))
			if(j==K) next
			for (i in (j+1):K){
				cl_pen <- cl_pen + fit.obj$gamma*sum(abs(fits[,j]-fits[,i])*fit.obj$subset.wts[[j]][[i-j]]*(sqrt(fit.obj$fit_var[,i] + fit.obj$fit_var[,j])))
			}
		}
		return(c(ssd, tf_pen, cl_pen))
}
