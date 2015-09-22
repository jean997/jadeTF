#Helper function for jade_path
find_new_gamma <- function(l1.total, log.gammas, sep.total, keep.fits,
                           start.step, l1.gap, l1.top, lg.top,
                           tol, buffer, verbose=TRUE){
  n <- length(keep.fits)
  #Too early to use smoothing
  if(n < 6){
    new.gamma <- max(log.gammas) + start.step
    return(new.gamma)
  }

  log.gammas.k <- log.gammas[keep.fits]
  l1.total.k <- l1.total[keep.fits]
  sep.total.k <- sep.total[keep.fits]
  #Top of the path
  if(l1.top == min(l1.total.k)){
    #We haven't really started into the part of the path we care about
    if(verbose) cat("Top!\n")
    l1.target <- min(l1.total.k) - l1.gap #Only one target
    new.gamma <- project_new_gamma_spline(l1.target=l1.target, lg.top=lg.top,
                                   l1.total=l1.total, log.gammas=log.gammas,
                                   l1.gap=l1.gap, buffer=buffer, keep=keep.fits)
    if(new.gamma$warn){
      #cat("Hmmm. Check what is happening\n")
      new.gamma$new.gamma <- max(log.gammas) + buffer
    }
    new.gamma <- new.gamma$new.gamma
    return(new.gamma)
  }

  lt <- sort(c(l1.total.k[ l1.total.k <= l1.top], 0), decreasing=TRUE)
  #Bottom of the path
  if(all(-1*diff(lt) <= (l1.gap*2))){
    #We have fits between l1.total0 and 0 with no gaps larger than 2*l1.gap. We might be done!
    if(verbose) cat("Bottom?\n")
    if( min(l1.total.k) < tol | min(sep.total.k)==0){
      #Close enough to the bottom. Done!
      return(NA)
    }else{
      #If the path is dense but we didn't get close to the bottom
      l1.target <- min(l1.total.k)/2 #Only one target
      new.gamma <- project_new_gamma_spline(l1.target=l1.target, lg.top=lg.top,
                                     l1.total=l1.total, log.gammas=log.gammas,
                                     l1.gap=l1.gap, buffer=buffer, keep=keep.fits)
      new.gamma <- new.gamma$new.gamma
      return(new.gamma)
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
    return(NA)
  }
  new.gamma <- new.gamma$new.gamma
  return(new.gamma)
}
