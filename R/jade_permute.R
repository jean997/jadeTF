#'Run JADE with permuted group assignments
#'@description This function takes the complete data matrix and re-runs
#'JADE permuting the group labels. It can optionally recalculate the weights
#'in the likelihood.
#'@param Y Full data matrix (p by n)
#'@param fit0 A JADE fit using original group labels at gamma=0
#'@param gammas List of gammas at which to fit permutation fits
#'@param tol Tolerance for declaring sites separated
#'@param n.perm Number of permutations
#'@param n.rep.fit.var If NULL use original fit.var. Otherwise
#'the number of bootstrap replications to use.
#'@param sd.type If "Orig" use sds used in fit0. Otherwsie recalculated sds
#'for each permutation assuming binomial or poisson models.
#'@param READS Required if sd.type="Binomial"
#'@return A list with elements
#'#' \describe{
#'  \item{\code{sep.total}}{A matrix length(gammas) by h n.perm }
#'  \item{\code{gammas}}{Gamma values fit for each permutation}
#'   \item{\code{tol}}{The tolerance at which the seaparation in \code{sep.total} and \code{sep}
#'   were calculated}
#' }
#'@export
jade_permute <- function(Y, fit0, gammas, save.prefix,
                         tol=1e-3, which.perm=1:100,
                         n.rep.fit.var=NULL, adjust.rho.alpha=TRUE,
                         sd.type=c("Orig", "Binomial", "Poisson"), READS=NULL){
  sd.type <- match.arg(sd.type)
  if(sd.type == "Binomial" & is.null(READS)){
    stop("For Binomial type please provide READS")
  }
  if(class(fit0)=="character"){
    fit0.file <- fit0
    fit0 <- getobj(fit0.file)
  }

  #For now, only admm supported here
  alg <- fit0$algorithm
  stopifnot(alg == "admm")

  #Info about the data
  n <- dim(Y)[2]
  p <- dim(Y)[1]
  sample.size <- fit0$sample.size
  K <- length(sample.size)
  grp.start <- c(1, cumsum(sample.size)[-K]+1)
  grp.stop <- cumsum(sample.size)
  if(!gammas[1]==0) gammas <- c(0, gammas)

  #Results to return
  sep.total <- matrix(nrow=length(gammas), ncol= length(which.perm))

  for(j in which.perm){
    #New data
    o <- sample(1:n, size = n, replace = FALSE)
    Y.perm <- Y[, o]
    if(!is.null(READS)) READS.perm <- READS[, o]
      else READS.perm <- NULL
    y.perm <- matrix(nrow=p, ncol=K)

    #Calc y.perm and new sds if nescessary
    sds.perm <- matrix(nrow=p, ncol=K)
    for(i in 1:K){
      Yi <- Y.perm[, grp.start[i]:grp.stop[i]]
      if(sd.type == "Orig" & is.null(READS)){
        y.perm[,i] <- rowSums(Yi)/sample.size[i]
      }else if(sd.type=="Binomial" | !is.null(READS)){
        Ri <- READS.perm[, grp.start[i]:grp.stop[i]]
        r <- binom.func(Yi, Ri)
        y.perm[,i] <- r$y
        if(sd.type=="Binomial") sds.perm[,i] <- r$sds
      }else if(sd.type=="Poisson"){
        r <- poisson.func(Yi)
        y.perm[,i] <- r$y
        sds.perm[,i] <- r$sds
      }
    }
    if(sd.type=="Orig") sds.perm <- fit0$sds
    #Fitvar
    if(!is.null(n.rep.fit.var)){
      fit.var <- matrix(nrow=p, ncol=K)
      for(i in 1:K){
        Yi <- Y.perm[, grp.start[i]:grp.stop[i]]
        if(!is.null(READS)) Ri <- READS.perm[, grp.start[i]:grp.stop[i]]
          else Ri <- NULL
        fit.var[,i] <- bs_var_tf(Y=Yi, lambda=fit0$lambda[i], sample.size=sample.size[i],
                            sd.type=sd.type, READS=Ri,
                            pos=fit0$pos, scale.pos=fit0$scale.pos, ord=fit0$ord,
                            n.rep=n.rep.fit.var)
      }
    }else{
      fit.var=fit0$fit.var
    }

    #Fit 0
    fit0.perm <- jade_admm(y=y.perm, gamma=0, pos=fit0$pos, scale.pos=fit0$scale.pos,
                     lambda=fit0$lambda, sample.size=sample.size, ord=fit0$ord,
                     sds=sds.perm, fit.var=fit.var)

    fn <- paste0(save.prefix, ".perm", j, ".RData")
    path.perm <- jade_path_planned(fit0.perm, out.file=fn,
                                  gammas=gammas[-1], return.object=TRUE,
                                  tol=tol, verbose = TRUE,
                                  adjust.rho.alpha=adjust.rho.alpha)
    s <- path.perm$sep.total
    if(!min(s) == 0) cat("Warning: Path may be incomplete.\n")
    my.idx <- match(round(path.perm$gammas, digits=8), round(gammas, digits=8))
    if(length(my.idx) != length(path.perm$gammas)) cat("Warning: Gammas may not match correctly")
    if(!all.equal(my.idx, 1:length(my.idx))) cat("Warning: Gammas may not match correctly")
    sep.total[my.idx, j] <- s
    sep.total[gammas > max(path.perm$gammas), j] <- 0
  }
  return(list("sep.total"=sep.total, "tol"=tol, "gammas"=gammas))
}

jade_permute_results <- function(save.prefix, orig.gammas, n.perm, new.tol=1e-3){
  i <- 1
  sep.total <- matrix(nrow=length(orig.gammas), ncol=n.perm)
  while(i <=n.perm){
    cat(i, " .. ")
    fn <- paste0(save.prefix, ".perm", i, ".RData")
    path.perm <- getobj(fn)
    if(new.tol){
      s <- get_sep_total(path.perm, tol=new.tol)
    }else{
      s <- path.perm$sep.total
    }
    if(!min(s) == 0) cat("Warning: Path may be incomplete.\n")
    my.idx <- match(path.perm$gammas, orig.gammas)
    if(length(my.idx) == length(path.perm$gammas)) cat("Warning: Gammas may not match correctly")
    sep.total[my.idx, i] <- s
    sep.total[orig.gammas > max(path.perm$gammas), i] <- 0
  }
  cat("\n")
  return(sep.total)
}


binom.func <- function(Y, READS){
  y <- rowSums(Y)/rowSums(READS)
  yhat <- (rowSums(Y)+0.5)/(rowSums(READS)+1)
  sds <- sqrt(yhat*(1-yhat)/rowSums(READS))
  return(list("y"=y, "sds"=sds))
}

poisson.func <- function(Y){
  n <- dim(Y)[2]
  y_total <- rowSums(Y)
  y <- y_total/n
  sds <- jitter(sqrt(y_total))/n
  return(list("y"=y, "sds"=sds))
}
