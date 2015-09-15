

jade_permute <- function(Y, fit0, gammas, save.prefix,
                         tol=1e-3, n.perm=100,
                         calc.fit.var =TRUE, n.rep.fit.var=100,
                         sd.type=c("Equal", "Binomial", "Poisson"), READS=NULL){
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

  n <- dim(Y)[2]
  p <- dim(Y)[1]
  sample.size <- fit0$sample.size
  K <- length(sample.size)
  grp.start <- c(1, cumsum(sample.size)[-K]+1)
  grp.stop <- cumsum(sample.size)

  i <- 1
  while(i <= n.perm){
    #New data
    order <- sample(1:n, size = n, replace = FALSE)
    Y.perm <- Y[, order]
    if(!is.null(READS)) READS.perm <- READS[, order]
      else READS.perm <- NULL
    y.perm <- matrix(nrow=p, ncol=K)
    sds.perm <- matrix(1, nrow=p, nocl=K)
    for(i in 1:K){
      Yi <- Y.perm[, grp.start[i]:grp.stop[i]]
      if(sd.type == "Equal"){
        y.perm[,i] <- rowSums(Yi)/sample.size[i]
      }else if(sd.type=="Binomial"){
        Ri <- READS.perm[, grp.start[i]:grp.stop[i]]
        r <- binom.func(Yi, Ri)
        y.perm[,i] <- r$y
        sds.perm[,i] <- r$sds
      }else if(sd.type=="Poisson"){
        r <- poisson.func(Yi)
        y.perm[,i] <- r$y
        sds.perm[,i] <- r$sds
      }
    }
    #Fitvar
    if(calc.fit.var){
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

    fn <- paste0(save.prefix, ".perm", i, ".RData")
    path.perm <- jade_path_planned(fit0.perm, out.file=fn,
                                  gammas=gammas, return.object=TRUE,
                                  tol=tol, verbose = TRUE,
                                  adjust.rho.alpha=TRUE)
  }
}

jade_permute_results <- function(save.prefix, orig.gammas, n.perm, new.tol=NULL){
  i <- 1
  sep.total <- matrix(nrow=length(orig.gammas), ncol=n.perm)
  while(i <=n.perm){
    cat(i, " .. ")
    fn <- paste0(save.prefix, ".perm", i, ".RData")
    path.perm <- getobj(fn)
    if(new.tol){
      s <- get_sep_total(path.perm, new.tol=new.tol)
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
