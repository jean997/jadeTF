#' Get a table of separated regions.
#'
#' @param fit A JADE fit.
#' @param which.window,chr Will be put in the table. Can be useful for writing other files.
#' @param new.tol Recalculate separation using a different tolerance.
#' @return A table of separated regions. Returns 0 if all profiles are fused or only
#' separated at singleton sites.
get_separated_regions <- function(fit, which.window=1, chr="chr22", new.tol=NULL, data.range=NULL){
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
  z.merge <- reduce(z, min.gapwidth=2)
  z.merge <- z.merge[ width(z.merge) > 1, ]

  #Keep partition types separated
  z.sep <- merge_close_regions(z)
  z.sep <- z.sep[ width(z.sep) > 1, ]

  if(nrow(z.merge) ==0 & nrow(z.sep) ==0) return(0)

  res.tab.sep <- data.frame(matrix(nrow=nrow(z.sep), ncol=7+(K*(K-1)/2)))
  names(res.tab.sep)[1:7] <- c("Chrom", "Start", "Stop", "Width", "Ncpg", "Window", "Partition")
  for(i in 8:ncol(res.tab.sep)){
    names(res.tab.sep)[i] <- paste("AvgGap", paste(idx_to_pair(i-7, K), collapse=""), sep="")
  }

  res.tab.merge <- data.frame(matrix(nrow=nrow(z.merge), ncol=6+(K*(K-1)/2)))
  names(res.tab.merge) <- c(names(res.tab.sep)[1:6], names(res.tab.sep)[8:(7+(K*(K-1)/2))])

  if(nrow(z.sep) > 0){
    res.tab.sep$Chrom <- chr
    res.tab.sep$Window <- which.window

    res.tab.sep$Start <- fit$pos[start(z.sep)]
    res.tab.sep$Stop <- fit$pos[end(z.sep)]
    res.tab.sep$Width <- res.tab.sep$Stop - res.tab.sep$Start
    res.tab.sep$Ncpg <- width(z.sep$ranges)

    for(i in 1:nrow(z.sep)){
      f <- fit$fits[start(z.merge)[i]:end(z.merge)[i],]
      p <- fit$pos[start(z.merge)[i]:end(z.merge)[i]]

      if(!is.null(data.range)){
        f <- pmin(f, max(data.range))
        f <- pmax(f, min(data.range))
      }
      res.tab.sep$Partition[i] <- sep_to_partition(z.sep$score[i], K)
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
    res.tab.merge$Width <- res.tab.merge$Stop - res.tab.merge$Start
    res.tab.merge$Ncpg <- width(z.merge$ranges)
    for(i in 1:nrow(z.merge)){
      for(j in 1:(K-1)){
        for(l in (j+1):K){
          g <- f[,j]-f[,l]
          res.tab.sep[i, paste("AvgGap", j, l, sep="")] <- sum(diff(p)*rollmean(abs(g),2))/(max(p)-min(p))
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
