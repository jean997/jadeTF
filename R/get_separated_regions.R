#' Get a table of separated regions.
#'
#' @param fit A JADE fit.
#' @param which.window,chr Will be put in the table. Can be useful for writing other files.
#' @param new.tol Recalculate separation using a different tolerance.
#' @param min.gapwidth Merge regions separated by fewer than this many sites
#' @param min.width Don't report regions smaller than this
#' @return A list of two tables of separated regions.
#' One of the tables is merged over partition types, the other is not.
#' Returns 0 if all profiles are fused or only. These tables ca be passed
#' to the \code{\link{plot_jade}} function in \code{sep.tab} argument.
#' separated at singleton sites.
#' @export
get_separated_regions <- function(fit, which.window=1, chr="chr22", new.tol=NULL, data.range=NULL,
                                  min.gapwidth=2, min.width=1){
  K <- dim(fit$fits)[2]
  #stopifnot(K==3)

  if(!is.null(new.tol)) sep <- get_sep(fit$beta, new.tol)
    else sep <- fit$sep

  sep <- matrix(unlist(sep), ncol=(K*(K-1)/2))

  if(all(sep==0)) return(0)

  S <- apply(sep, MARGIN = 1, FUN = paste, collapse=" ")
  zeros <- paste(rep(0, (K*(K-1)/2)), collapse=" ")

  z <- as(Rle(S), "RangedData")
  if(all( z$score == zeros)) return(0)

  #Get rid of fully fused regions
  z <- z[ !score(z) == zeros, ]

  #Merge over partition types
  z.merge <- reduce(z, min.gapwidth=min.gapwidth)
  z.merge <- z.merge[ width(z.merge) > min.width, ]

  #Keep partition types separated
  z.sep <- merge_close_regions(z, margin=min.gapwidth)
  z.sep <- z.sep[ width(z.sep) > 1, ]

  if(nrow(z.merge) ==0 & nrow(z.sep) ==0) return(0)

  res.tab.sep <- data.frame(matrix(nrow=nrow(z.sep), ncol=6+(K*(K-1)/2)))
  names(res.tab.sep)[1:7] <- c("Chrom", "Start", "Stop", "Width", "Window", "Partition")
  for(i in 8:ncol(res.tab.sep)){
    names(res.tab.sep)[i] <- paste("AvgGap", paste(idx_to_pair(i-6, K), collapse=""), sep="")
  }

  res.tab.merge <- data.frame(matrix(nrow=nrow(z.merge), ncol=5+(K*(K-1)/2)))
  names(res.tab.merge) <- c(names(res.tab.sep)[1:5], names(res.tab.sep)[8:(6+(K*(K-1)/2))])

  if(nrow(z.sep) > 0){
    res.tab.sep$Chrom <- chr
    res.tab.sep$Window <- which.window

    res.tab.sep$Start <- fit$pos[start(z.sep)]
    res.tab.sep$Stop <- fit$pos[end(z.sep)]
    res.tab.sep$Width <- res.tab.sep$Stop - res.tab.sep$Start + 1

    for(i in 1:nrow(z.sep)){
      f <- fit$fits[start(z.sep)[i]:end(z.sep)[i],]
      p <- fit$pos[start(z.sep)[i]:end(z.sep)[i]]

      if(!is.null(data.range)){
        f <- pmin(f, max(data.range))
        f <- pmax(f, min(data.range))
      }
      my.sep <- as.numeric(unlist(strsplit(z.sep$score[i], split=" ")))
      res.tab.sep$Partition[i] <- sep_to_partition(my.sep, K, data=f)
      for(j in 1:(K-1)){
        for(l in (j+1):K){
          g <- f[,j]-f[,l]
          res.tab.sep[i, paste("AvgGap", j, l, sep="")] <- sum(diff(p)*rollmean(abs(g),2))/(max(p)-min(p))
        }
      }
    }
  }

  if(nrow(z.merge) > 0){
    res.tab.merge$Chrom <- chr
    res.tab.merge$Window <- which.window
    res.tab.merge$Start <- fit$pos[start(z.merge)]
    res.tab.merge$Stop <- fit$pos[end(z.merge)]
    res.tab.merge$Width <- res.tab.merge$Stop - res.tab.merge$Start + 1
    for(i in 1:nrow(z.merge)){
      f <- fit$fits[start(z.merge)[i]:end(z.merge)[i],]
      p <- fit$pos[start(z.merge)[i]:end(z.merge)[i]]

      if(!is.null(data.range)){
        f <- pmin(f, max(data.range))
        f <- pmax(f, min(data.range))
      }
      for(j in 1:(K-1)){
        for(l in (j+1):K){
          g <- f[,j]-f[,l]
          res.tab.merge[i, paste("AvgGap", j, l, sep="")] <- sum(diff(p)*rollmean(abs(g),2))/(max(p)-min(p))
        }
      }
    }
  }
  return(list("merged"=res.tab.merge, "separated"=res.tab.sep))

}

#Merge preserving partition type
#z should be a RangedData object with a value collumn called "score"
merge_close_regions <- function(z, margin=2){
  stopifnot(class(z)=="RangedData")
  starts <- c()
  ends <- c()
  scores <- c()
  ss <- unique(score(z))
  for(s in ss){
    q <- z[score(z) == s,]
    q <- reduce(q, min.gapwidth=margin)
    ends <- c(ends, end(q))
    starts <- c(starts, start(q))
    scores <- c(scores, rep(s, nrow(q)))
  }
  j <- order(starts)
  ir <- IRanges(starts[j], ends[j])
  rd <- RangedData(ir)
  rd$score <- scores[j]
  return(rd)
}
