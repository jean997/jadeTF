#Methylation data utility functions

#' Combine data across replicates
#'
#' @param dfs A list of data frames. Each data frame
#' should have three columns - position, reads, and counts
#' @return A combined data frame with three columns - position (pos), number of reads (reads), and
#' number of counts (y)
#' @export
combine_replicates <- function(dfs, range=NULL){
  #Add together replicates
  #dfs have collumns pos, reads, y
  P <- c() #Positions in any of the dfs
  n <- length(dfs)
  cat(n, "\n")
  #get positions
  for(j in 1:n){
    pos <- dfs[[j]]$pos
    P <- sort(unique(c(P, pos)))
  }
  cat("Positions: ", length(P), "\n")
  if(!is.null(range)){
    P <- P[P <= max(range) & P >= min(range)]
    cat("Positions remaining: ", length(P), "\n")
  }

  y <- rep(0, length(P))
  reads <- rep(0, length(P))
  for(j in 1:n){
    cat(j, "\n")
    ind <- which(P %in% dfs[[j]]$pos)
    data.ind <- which(dfs[[j]]$pos %in% P)
    cat(length(ind), " " , length(data.ind), "\n")
    y[ind] <- y[ind] + dfs[[j]]$counts[data.ind]
    reads[ind] <- reads[ind] + dfs[[j]]$reads[data.ind]
  }
  return(data.frame("counts"=y, "reads"=reads, "pos"=P))
}

#' Take multiple data frames with data from a single condition and
#'make data ready for JADE
#'
#' @param dfs A list of data frames. Each data frame
#' should have three columns - position, reads, and counts
#' @return A list with elements counts, reads, pos, probs, and sds.
#' @export
collect_methylation_data <- function(dfs){
  #Combine multiple data sets, format for jade
  K <- length(dfs)
  P <- c()
  #Get positions
  for(j in 1:K){
    pos <- dfs[[j]]$pos
    P <- sort(unique(c(P, pos)))
  }
  p <- length(P)
  y <- matrix(NA, p, K)
  reads <- matrix(NA, p, K)
  probs <- matrix(NA, p, K)
  sds <- matrix(NA, p, K)
  for (j in 1:K){
    idx <- which(P %in% dfs[[j]]$pos)
    y[idx,j] <- dfs[[j]]$counts
    reads[idx,j] <- dfs[[j]]$reads
    probs[idx,j] <- y[idx,j]/reads[idx,j]
    phat <- (y[idx,j] + 0.5)/(reads[idx,j]+1)
    sds[idx,j] <- sqrt(phat*(1-phat)/reads[idx,j])
  }
  return(list("counts"=y, "reads"=reads, "pos"=P, "probs"=probs, "sds"=sds))
}

#' Find window boundaries for methylation data
#'
#' @param positions Vector of positions
#' @param gapsize Maximum gap to allow between data points in a window
#' @param mindata Minimum number of data points per window
#' @param nm.edges Make sure the last data point in each window is non-missing in all samples
#' @param nmiss If nm.edges=TRUE provide a vector listing the number
#' of missing points at each position.
#' @return A list with elements counts, reads, pos, probs, and sds.
#' @export
methylation_windows <- function(positions, gapsize=2000,
                                mindata=20, nmiss=NULL, nm.edges=FALSE){
  p <- length(positions)
  diff <- positions[2:p]-positions[1:(p-1)]
  z <- rle(diff > gapsize)
  nsegs <- sum(z$lengths > mindata & z$values==FALSE)
  K <- which(z$lengths > mindata & z$values ==FALSE)
  ncpg <- z$lengths[K]+1
  start <- c()
  stop <- c()
  #cat(length(K), "\n")
  for(k in K){
    if(k==1){ start <- c(start, 1);}
    else{ start<- c(start, sum(z$lengths[1:(k-1)])+1)}
    stop <- c(stop, sum(z$lengths[1:k])+1)
  }

  if(nm.edges){
    if(is.null(nmiss)) stop("Provide missing counts")
    for(i in 1:nsegs){
      nm <- nmiss[start[i]:stop[i]]
      l <- length(nm)
      if(all(nm > 0)){
        cat("Remove ", i, "\n")
      }else if(nm[1] ==0 & nm[l] ==0){
        next
      }
      fct <- 0; bct <- 0
      while(nm[1] > 0){
        start[i] <- start[i] +1
        ncpg[i] <- ncpg[i]-1
        nm <- nmiss[start[i]:stop[i]]
        l <- length(nm)
        fct <- fct+1
      }
      while(nm[l] > 0){
        stop[i] <- stop[i]-1
        ncpg[i] <- ncpg[i]-1
        nm <- nmiss[start[i]:stop[i]]
        l <- length(nm)
        bct <- bct+1
      }
      #cat(i, ": ", fct, " ", bct, "\n")
    }
  }
  start.pos <- positions[start]; stop.pos <- positions[stop]
  size <- stop.pos-start.pos
  df <- data.frame(start, stop, ncpg, start.pos, stop.pos, size)
  df <- df[ncpg >= mindata,]
  return(df)
}
