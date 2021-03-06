#' @useDynLib jadeTF
#' @importFrom Rcpp sourceCpp
#' @import  genlasso
#' @import  clusterpathRcpp
#' @import Matrix
#' @import zoo
#' @import plotrix
#' @import IRanges
#' @import ggplot2

#Contains small utility functions

#Default weights
default_wts <- function(p, K){
  wts <- list()
  if(K==1){
    wts[[1]] <- rep(1, p)
    return(wts)
  }
  for(j in 1:(K-1)){
    wts[[j]] <- list()
    for(i in (j+1):K){
      wts[[j]][[i-j]] <- rep(1, p)
    }
	}
	return(wts)
}

wts_from_var <- function(fit.var){
  wts <- list()
  for(j in 1:(K-1)){
    wts[[j]] <- list()
    for(i in (j+1):K){
      wts[[j]][[i-j]] <- sqrt(fit.var[,i] + fit.var[,j])
    }
  }
  return(wts)
}

pairwise_wts <- function(subset.wts, fit.var, sample.size){
	K <- length(subset.wts)+1
	pw <- list()
	for(j in 1:(K-1)){
		pw[[j]] <- list()
		for(i in (j+1):K){
			pw[[j]][[i-j]] <- subset.wts[[j]][[i-j]]*sqrt( fit.var[i] + fit.var[j])
		}
	}
	return(pw)
}


#For jade_gd: Duals get initialized on the boundary if initial values not given
starting_duals <- function(fit, gamma, var.wts, sample.size, subset.wts){
  p <- dim(fit)[1]
  K <- dim(fit)[2]
  duals <- list()
  for(j in 1:(K-1)){
    duals[[j]] <- list()
    for(i in (j+1):K){
      u<- gamma*var.wts[[j]][[i-j]] * as.vector(sign(fit[,j]-fit[,i]))
      w <- subset.wts[[j]][[i-j]]
      u[w==0] <- 0
      duals[[j]][[i-j]] <- u
    }
  }
  return(duals)
}
constrain_duals <- function(duals, gamma, var.wts, sample.size, subset.wts){
  p <- length(duals[[1]][[1]])
  K <- length(sample.size)
  for(j in 1:(K-1)){
    for(i in (j+1):K){
      u<- duals[[j]][[i-j]]
      u <- sign(u)*pmin(abs(u), gamma*var.wts[[j]][[i-j]])
      w <- subset.wts[[j]][[i-j]]
      u[w==0] <- 0
      duals[[j]][[i-j]] <- u
    }
  }
  return(duals)
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
#' Get an object from an .RData file
#' @export
getobj <- function (Rdata){
    objname <- load(Rdata)
    if (length(objname) > 1) {
        warning(paste("Multiple objects stored in file", Rdata,
            "\nReturning only the first object"))
    }
    return(get(objname))
}


#' Determine which pairs of sites are separated.
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

#' Recalculate sep.total for a path that has already been run
#' @param path Path of JADE fits
#' @param tol Tolerance for determining separation
#' @return A vector of same length as path$gammas
#' @export
get_sep_total <- function(path, tol){
  sep.total <- rep(NA, length(path$JADE_fits))
  sep.total[1] <- sum( unlist(get_sep(path$JADE_fits[[1]]$fits, tol=tol)))
  s <- unlist( lapply(path$JADE_fits[-1],  FUN = function(f){
    sep <- sum(unlist(get_sep(f$beta, tol=tol)))
    return(sep)}))
  sep.total[-1] <- s
  return(sep.total)
}

#' Calculate total number of separated regions for each fit in a path
#' @param path Path of JADE fits
#' @param tol Tolerance for determining separation
#' @return A vector of same length as path$gammas
#' @export
get_regions_total <- function(path, tol){
  reg.total <- rep(NA, length(path$JADE_fits))
  sep1 <- rowSums(matrix(unlist(get_sep(path$JADE_fits[[1]]$fits, tol=tol))))
  reg.total[1] <- sum(rle(sep1)$values)
  r <- unlist( lapply(path$JADE_fits[-1],  FUN = function(f){
    sep <- rowSums(matrix(unlist(get_sep(f$beta, tol=tol))))
    reg <- sum(rle(sep)$values)
    return(reg)}))
  reg.total[-1] <- r
  return(reg.total)
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

plot_separation <- function(path.obj, sites, hline=NULL, main="", log=TRUE, ylim=NULL){
  K <- length(path.obj$sep)+1
  y =log10(path.obj$gammas[-1])
  k <- length(y)
  p <- length(sites)
  if(is.null(ylim)) ylim <- range(y)
  idx <- order(y, decreasing = FALSE)
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      sep <- path.obj$sep[[i]][[j-i]]
      image(x=sites, y=y[idx], z=sep[, idx], col=c("white", "blue"), main=main, ylim=ylim)
      if(!is.null(hline)) abline(h=hline, col="red")
    }
  }
}

#Objective

#JADE Objective function for trend filtering problem
#These are used by jade_gd
obj_fct_wrapper <- function(jade.obj){

  if(!is.null(jade.obj$scale.pos)){
    pos <- jade.obj$pos
    R <-range(pos)
    pos<- jade.obj$scale.pos* ((pos-R[1])/(R[2]-R[1]))
    jade.obj$pos <- pos
  }

  obj.value <- obj_fct(jade.obj$y, jade.obj$fits, jade.obj$lambda, jade.obj$gamma, jade.obj$sample.size,
           jade.obj$sds, jade.obj$pos, jade.obj$ord)
  return(obj.value)
}
obj_fct <-  function(y, theta, lambda, gamma, sample.size,
                     sds,  pos, ord){

  h <- function(x, p, pos, ord){
    D_tf <- getDtfPosSparse(n=p, ord=ord, pos=pos)
    return(sum(abs ( D_tf %*% x)))
  }

  p <- dim(y)[1]
  K <- dim(y)[2]

  obj.value <- 0

  for(j in 1:K){
    nm <- which(!is.na(y[,j]))
    obj.value <- obj.value + (sample.size[j]/2)*sum((1/sds[nm,j])*(y[nm,j]-theta[nm,j])^2) +
      lambda[j]*h(theta[,j], p, pos, ord)
    if(j==K) next
    for (i in (j+1):K){
      pen <- gamma*sum(abs(theta[,j]-theta[,i]))
      obj.value <- obj.value + pen
    }
  }
  return(obj.value)
}

expit <- function(x){return(exp(x)/(1 + exp(x)))}
logit <- function(x){return( log(x/(1-x)))}
