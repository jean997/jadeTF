#Helper function for jade_path
find_new_gamma <- function(l1.total, log.gammas, sep.total, n.fits,
                           converged, start.step, tol, buffer, verbose=TRUE){
  sep.total0 <- sep.total[1]
  if(all(sep.total[-1] < sep.total0)) sep.total0 <- sep.total[2]
  l1.total0 <- l1.total[1]

  #Only base next gamma on fits with l1.total < l1.total0 (unless < 6)
  keep.fits <- which(l1.total <= l1.total0 & is.finite(log.gammas) & converged)
  #Too early to use smoothing
  if(length(keep.fits) < 6){
    #Make start.step larger if many steps and no convergence
    f <- floor((sum(l1.total > l1.total0 | !converged))/3)
    start.step <- start.step*(2^f)
    new.gamma <- max(log.gammas) + start.step
    return(list("new.gamma"=new.gamma, "l1.gap"=min(l1.total)/n.fits, "done"=FALSE))
  }

  log.gammas.k <- log.gammas[keep.fits]
  l1.total.k <- l1.total[keep.fits]
  sep.total.k <- sep.total[keep.fits]

  if(any(sep.total.k >= sep.total0)){
    l1.top <- min(l1.total.k[sep.total.k >= sep.total0])
    lg.top <- max(log.gammas.k[sep.total.k >= sep.total0])
  }else{
    l1.top <- max(l1.total.k)
    lg.top <- min(log.gammas.k)
  }
  l1.gap <- l1.top/n.fits

  #Top of the path
  if(l1.top == min(l1.total.k)){
    #We haven't really started into the part of the path we care about
    if(verbose) cat("Top!\n")
    l1.target <- min(l1.total.k) - l1.gap #Only one target
    new.gamma <- project_new_gamma_spline(l1.target=l1.target, lg.top=lg.top,
                                   l1.total=l1.total, log.gammas=log.gammas,
                                   l1.gap=l1.gap, buffer=buffer, keep.fits=keep.fits)
    if(new.gamma$warn){
      #cat("Hmmm. Check what is happening\n")
      new.gamma$new.gamma <- max(log.gammas) + buffer
    }
    new.gamma <- new.gamma$new.gamma
    return(list("new.gamma"=new.gamma, "l1.gap"=l1.gap, "done"=FALSE))
  }

  lt <- sort(c(l1.total.k[ l1.total.k <= l1.top], 0), decreasing=TRUE)
  #Bottom of the path
  if(all(-1*diff(lt) <= (l1.gap*2))){
    #We have fits between l1.total0 and 0 with no gaps larger than 2*l1.gap. We might be done!
    if(verbose) cat("Bottom?\n")
    if( min(l1.total.k) < tol | min(sep.total.k)==0){
      #Close enough to the bottom. Done!
      return(list("new.gamma"=NA, "l1.gap"=l1.gap, "done"=TRUE))
    }else{
      #If the path is dense but we didn't get close to the bottom
      l1.target <- min(l1.total.k)/2 #Only one target
      new.gamma <- project_new_gamma_spline(l1.target=l1.target, lg.top=lg.top,
                                     l1.total=l1.total, log.gammas=log.gammas,
                                     l1.gap=l1.gap, buffer=buffer, keep=keep.fits)
      new.gamma <- new.gamma$new.gamma
      return(list("new.gamma"=new.gamma, "l1.gap"=l1.gap, "done"=FALSE))
    }
  }

  if(verbose) cat("Middle\n")
  #We are in the middle but the path isn't dense enough
  hole.idx <- which(-1*diff(lt) > 2*l1.gap)
  l1.target <- lt[hole.idx]-l1.gap #Possibly multiple targets
  if(min(l1.total.k) > tol | min(sep.total.k) > 0){
    if (all(l1.target > min(l1.total.k))){
      l1.target <- c(l1.target, min(l1.total.k)/2)
    }
  }
  new.gamma <- project_new_gamma_spline(l1.target=l1.target, lg.top=lg.top,
                                 l1.total=l1.total, log.gammas=log.gammas,
                                 l1.gap=l1.gap, buffer=buffer, keep=keep.fits)
  if(new.gamma$warn & (min(l1.total.k) < tol | min(sep.total.k) == 0)){
    cat("Warning: Path may have holes!\n")
    return(list("new.gamma"=NA, "l1.gap"=l1.gap, "done"=TRUE))
  }
  new.gamma <- new.gamma$new.gamma
  return(list("new.gamma"=new.gamma, "l1.gap"=l1.gap, "done"=FALSE))
}
