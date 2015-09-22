#Helper function for jade_path
#Function for finding the next gamma value give
#l1.target List of desired values of l1.total at the next fit
#lg.top
#l1.total List of l1.total values for previous fits
#log.gammas Corresponding list of log.gamma
#buffer
#l1.gap
project_new_gamma_spline <- function(l1.target, l1.total, log.gammas,
                              keep, lg.top, buffer, l1.gap){

  j <- length(keep)
  l1.total.k <- l1.total[keep]
  log.gammas.k <- log.gammas[keep]
  y=l1.total.k[order(log.gammas.k)]
  x=sort(log.gammas.k)

  #cat(length(x), ord, k, "\n")
  #cat(y, "\n", x, "\n")
  x.new <- seq(lg.top-2*buffer, max(log.gammas)+1, by=(buffer/2))
  sm <-try(smooth.spline(x=x, y=y), silent=TRUE)
  if(class(sm)=="try-error"){
    sm <- lm(y~x)
    pred <- approx(x = x, y=sm$fitted.values, xout = x.new)
  }else{
    pred <- predict(sm, x = x.new)
  }
  co.allgammas <- cummin(pred$y)
  new.gamma <- NA
  i <- 1
  warn <- FALSE
  while(is.na(new.gamma) & i <= j){
    new.gamma <- x.new[ which.min(abs(co.allgammas-l1.target[i]))]
    #Not allowed to go back too far
    new.gamma <- max(new.gamma, lg.top-2*buffer)

    #Don't repeat a gamma we have already tried
    ng <- new.gamma
    while(min(abs(log.gammas-new.gamma)) <= buffer) new.gamma <- new.gamma+buffer
    cat(l1.target[i], " ", new.gamma, "\n")
    if(abs(new.gamma-ng) > 5*buffer){
      #We haven't been able to hit this target - try the next one
      if(i < j) new.gamma <- NA
        else warn <- TRUE
    }
    i <- i+1
  }
  return(list("new.gamma"=new.gamma, "warn"=warn))
}
