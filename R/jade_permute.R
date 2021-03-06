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
jade_permute <- function(Y, fit0, gammas, save.prefix, sample.size=NULL,
                         tol=1e-3, which.perm=1:100, sds.mat=NULL,
                         n.rep.fit.var=NULL, adjust.rho.alpha=TRUE,
                         sd.type=c("Mat2", "Binomial", "Poisson", "EqualFromData", "Mat", "Orig"), READS=NULL){
  sd.type <- match.arg(sd.type)
  if(sd.type=="Mat" | sd.type=="Mat2"){
    if(is.null(sds.mat)) stop("For Mat sd.type please provide matrix sds.mat")
    stopifnot(all(dim(Y)==dim(sds.mat)))
    stopifnot(all(sds.mat > 0))
    V <- sds.mat^2
  }
  if(sd.type == "Binomial" & is.null(READS)){
    stop("For Binomial type please provide READS")
  }
  if(class(fit0)=="character"){
    fit0.file <- fit0
    fit0 <- getobj(fit0.file)
  }
  if(!sd.type=="Orig" & !all(fit0$sample.size == 1)) stop("Error: With sd.type ", sd.type, " sample size should be 1 for all groups.")
  if(is.null(sample.size)) sample.size <- fit0$sample.size
  stopifnot(all(gammas > 0))
  #For now, only admm supported here
  alg <- fit0$algorithm
  stopifnot(alg == "admm")

  #Info about the data
  n <- dim(Y)[2]
  p <- dim(Y)[1]
  stopifnot(sum(sample.size) == n)
  K <- length(sample.size)
  grp.start <- c(1, cumsum(sample.size)[-K]+1)
  grp.stop <- cumsum(sample.size)


  #Results to return
  sep.total <- matrix(nrow=length(gammas)+1, ncol= length(which.perm))

  for(j in 1:length(which.perm)){
    #New data
    o <- sample(1:n, size = n, replace = FALSE)
    Y.perm <- Y[, o]
    if(!is.null(sds.mat)){
      sds.mat.perm <- sds.mat[,o]
      V.perm <- V[,o]
    }
    if(!is.null(READS)) READS.perm <- READS[, o]
      else READS.perm <- NULL
    y.perm <- matrix(nrow=p, ncol=K)

    #Calc y.perm and new sds if nescessary
    sds.perm <- matrix(nrow=p, ncol=K)
    for(i in 1:K){
      Yi <- Y.perm[, grp.start[i]:grp.stop[i]]
      if(sd.type=="Mat"){
        sds.mati <- sds.mat.perm[, grp.start[i]:grp.stop[i]]
        r <- mat.func(Yi, sds.mati)
        y.perm[,i] <- r$y
        sds.perm[,i]<- r$sds
      }else if(sd.type=="Mat2"){
        Vi <- V.perm[, grp.start[i]:grp.stop[i]]
        r <- mat.func2(Yi, Vi)
        y.perm[,i] <- r$y
        sds.perm[,i]<- r$sds
      }else if(sd.type=="Orig"){
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
      }else if(sd.type=="EqualFromData"){
        r <- eqfd.func(Yi)
        y.perm[,i] <- r$y
        sds.perm[,i]<- r$sds
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
                            pos=fit0$pos, ord=fit0$ord,
                            n.rep=n.rep.fit.var)
      }
    }else{
      fit.var=fit0$fit.var
    }

    #Fit 0
    fit0.perm <- jade_admm(y=y.perm, gamma=0, pos=fit0$pos,
                     lambda=fit0$lambda, sample.size=fit0$sample.size, ord=fit0$ord,
                     sds=sds.perm, fit.var=fit.var)

    fn <- paste0(save.prefix, ".perm", which.perm[j], ".RData")

    path.perm <- jade_path_planned(fit0.perm, out.file=fn,
                                  gammas=gammas, return.object=TRUE,
                                  tol=tol, verbose = TRUE,
                                  adjust.rho.alpha=adjust.rho.alpha)
    s <- path.perm$sep.total
    my.idx <- match(round(path.perm$gammas, digits=8), round(c(0, gammas), digits=8))
    if(length(my.idx) != length(path.perm$gammas)) cat("Warning: Gammas may not match correctly")
    if(!all.equal(my.idx, 1:length(my.idx))) cat("Warning: Gammas may not match correctly")
    sep.total[ my.idx, j] <- s
    sep.total[c(0, gammas) > max(path.perm$gammas), j] <- 0
  }
  return(list("sep.total"=sep.total, "tol"=tol, "gammas"=gammas))
}

jade_permute_results <- function(save.prefix, gammas, n.perm, tol){
  i <- 1
  stopifnot(all(gammas > 0))
  sep.total <- matrix(nrow=length(gammas)+1, ncol=n.perm)
  while(i <=n.perm){
    cat(i, " .. ")
    fn <- paste0(save.prefix, ".perm", i, ".RData")
    path.perm <- getobj(fn)
    s <- get_sep_total(path.perm, tol=tol)

    my.idx <- match(round(path.perm$gammas, digits=8), round(c(0, gammas), digits=8))
    if(length(my.idx) != length(path.perm$gammas)) cat("Warning: Gammas may not match correctly")
    if(!all.equal(my.idx, 1:length(my.idx))) cat("Warning: Gammas may not match correctly")
    sep.total[ my.idx, i] <- s
    sep.total[c(0, gammas) > max(path.perm$gammas), i] <- 0
    i <- i+1
  }
  cat("\n")
  return(list("sep.total"=sep.total, "tol"=tol, "gammas"=gammas))
}

mat.func <- function(Y, sds.mat){
  S <- 1/(sds.mat^2)
  y <- rowSums(Y/(sds.mat^2))/rowSums(S)
  sds <- 1/sqrt(rowSums(S))
  return(list("y"=y, "sds"=sds))
}
mat.func2 <- function(Y, V){
  y <- rowMeans(Y)
  sds <- sqrt(rowSums(V))/ncol(V)
  return(list("y"=y, "sds"=sds))
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
  sds <- sqrt(y_total)/n
  return(list("y"=y, "sds"=sds))
}
eqfd.func <- function(Y){
  y <- rowMeans(Y)
  sds <- apply(Y, MARGIN=1, FUN=sd)/ncol(Y)
  sds[sds==0] <- min(sds[sds > 0])
  return(list("y"=y, "sds"=sds))
}

huber.func <- function(Y, k, min.var){
  yv <- apply(Y, MARGIN=1, FUN=function(x){
    mu <-  hubers(x, k=k, s=1)
    v <- max(huber_var(x, k=k, muH=mu), min.var)
    return(c(mu, v))
  })
  return(list("y"=yv[1,], "sds"=sqrt(yv[2,])))
}

huber_var <- function(x, muH, k){
  n <- length(x)
  ifs <- x-muH
  ifs[(x-muH) < -k] <- -k
  ifs[(x-muH) > k] <- k
  C <- sum(abs(x-muH)<=k)/n
  ifs <- ifs/C
  return(sum(ifs^2)/n^2)
}
