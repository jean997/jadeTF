#Helper function for jade_path
#Function for finding the next gamma value give
#l1.target List of desired values of l1.total at the next fit
#lg.top
#l1.total List of l1.total values for previous fits
#log.gammas Corresponding list of log.gamma
#buffer
#l1.gap
project_new_gamma_smooth <- function(l1.target, l1.total, log.gammas,
                              lg.top, buffer, l1.gap){

  j <- length(l1.target)

  y=l1.total[order(log.gammas)]
  x=sort(log.gammas)
  if(length(x) < 15) ord <- 1
  else ord <- 2

  if(length(x) < 25) k <- 3
  else k <- 5

  #cat(length(x), ord, k, "\n")
  #cat(y, "\n", x, "\n")
  sm <- genlasso::trendfilter(pos=x, y=y, ord=ord)
  sm.cv <- cv.trendfilter(sm, k=k)
  co <- coef.genlasso(sm, lambda=sm.cv$lambda.1se)
  x.new <- seq(lg.top-2*buffer, max(log.gammas)+1, by=(buffer/2))
  co.allgammas = .Call("tf_predict_R",
             sBeta = as.double(co$beta),
             sX = as.double(x),
             sN = length(y),
             sK = as.integer(ord),
             sX0 = as.double(x.new),
             sN0 = length(x.new),
             sNLambda = 1,
             sFamily = 0,
             sZeroTol = as.double(1e-6),
             package = "jadeTF")

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

