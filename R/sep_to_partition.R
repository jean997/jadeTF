#Utility functions for dealing with partitions

sep_to_partition <- function(sep, K, data=NULL){
	if(all(sep==0)){
		partition <- paste(1:K, collapse=" ")
		partition <- paste("(", partition, ")", sep="")
		return(partition)
	}
	if(all(sep==1)){
	  if(!is.null(data)){
	    avg.value <- colMeans(data)
	    o <- order(avg.value)
	    partition <- paste((1:K)[o], collapse=")(")
	  }else{
	     partition <- paste(1:K, collapse=")(")
	  }
		partition <- paste("(", partition, ")", sep="")
		return(partition)
	}

  pairs <- unlist(lapply(which(sep==0), FUN=function(x){jadeTF:::idx_to_pair(x, K)}))
	pairs <- matrix(pairs, ncol=2, byrow = TRUE)


	groups <- list(c(1))
	for(i in 1:nrow(pairs)){
		p <- pairs[i,]
		found=FALSE
		for(j in 1:length(groups)){
			if(any(p %in% groups[[j]])){
				groups[[j]] <- unique(c(groups[[j]], p))
				found <- TRUE
				break
			}
		}
		j <- length(groups)
		if(!found) groups[[j+1]] <- p
	}
  for(i in 1:K){
    found <- FALSE
    for(j in 1:length(groups)){
      if(i %in% groups[[j]]) found <- TRUE
    }
    if(!found) groups[[j+1]] <- c(i)
  }

  #Sort groups
  if(!is.null(data)){
    avg.value <- rep(NA, length(groups))
    for(j in 1:length(groups)){
      p <- groups[[j]][1]
      avg.value[j] <- mean(data[,p], na.rm=TRUE)
    }
    o <- order(avg.value)
    groups.sorted <- list()
    for(j in 1:length(groups)){
      groups.sorted[[j]] <- groups[[o[j]]]
    }
    groups <- groups.sorted
  }
	partition <- ""
	for(j in 1:length(groups)){
		pt <- paste("(", paste( groups[[j]], collapse=" "), ")", sep="")
		partition <- paste(partition, pt, sep="")
	}

	exp.sep <- partition_to_sep(partition, K)
	if(!all(sep == exp.sep)){
		cat("WARNING: Inconsistent sep vector\n")
		cat("Expected: ", exp.sep, "\nGiven", sep, "\n")
		partition <- paste0("**", partition)
	}
	return(partition)
}


partition_to_sep <- function(partition, K){
	grps <-unlist(strsplit(partition, split="[()]"))
	grps <- grps[!grps==""]
	sep <- rep(1, K*(K-1)/2)
	for(j in 1:length(grps)){
		grp <- as.numeric(unlist(strsplit(grps[j], split=" ")))
    if(length(grp) == 1) next
		for(i in 1:(length(grp)-1)){
			for(k in (i+1):length(grp)){
				#cat(grp[i], grp[k], "\n")
				idx <- pair_to_idx(grp[i], grp[k], K)
				sep[idx] <- 0
			}
		}
	}
	return(sep)
}
