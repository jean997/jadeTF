#' @useDynLib jadeTF
#' @importFrom Rcpp sourceCpp
#' @import  genlasso
#' @import  clusterpathRcpp
#' @import  Matrix
#' @import zoo
#' @import plotrix
#' @import IRanges

#Contains small utility functions

#Fit one sample
#Missing points are imputed
#Minimize N/(2) || (y - \theta)/sds ||^2 + \lambda_1||D\theta||_1
#Equivalent to 1/2 || w*(y - \theta) ||^2 + \lambda_1/N||D\theta||_1
fit_one <- function(y, lambda, pos, sds, sample.size, ord){
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
      cv <- cv_pred.genlasso(obj=tfit.out, n.folds = 5, mode = "predict")
      l <- cv$lambda.1se
      lambda <- l*sample.size
      if(equal.wts) lambda <- lambda*(wts[1])^2
      cat(lambda, "\n")
		}else if(equal.wts){
			l <- lambda/(sample.size*wts[1]^2)
		}else{
			l <- lambda/sample.size
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
						lambda[j] <- f$lambda
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

		fit <- fit_one(new.y, new.lam, pos, new.sigma, 1, ord)
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


#' Determine which pairs of sites are fused.
#' @param fits A p x K matrix of fits. To recalculate separation
#' for a \code{jade_admm} object use \code{fits=obj$beta}.
#' @param tol Tolerance for determining separation
#' @return A list of lists of length p vectors containing 0s and 1s.
#' The vector stored in \code{[[i]][[j-i]]} indicates the separation between group i and group j.
#' @export
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

pair_to_idx <- function(i, j, K){
  if(i==j) stop("i == j given to pair_to_idx")
  my.i <- min(i, j)
  my.j <- max(i, j)

  idx <- 1
  first <- 1
  while(my.i > first){
    idx <- idx + (K-first)
    first <- first + 1
  }
  idx <- idx + (my.j - my.i)-1
  return(idx)
}

idx_to_pair <- function(idx, K){
  my.i <- 1
  my.idx <- idx
  while(my.idx > 0){
    my.i <- my.i + 1
    my.idx <- my.idx - (K-my.i)-1
  }
  my.j <- my.idx + (K-my.i) +1 + (my.i -1)
  my.i <- my.i -1
  return(c(my.i, my.j))
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


#Objective

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

expit <- function(x){return(exp(x)/(1 + exp(x)))}
logit <- function(x){return( log(x/(1-x)))}
