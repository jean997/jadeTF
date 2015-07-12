
sep_to_partition <- function(sep, K){
	if(all(sep==0)){
		partition <- paste(1:K, collapse=" ")
		partition <- paste("(", partition, ")", sep="")
		return(partition)
	}
	if(all(sep==1)){
		partition <- paste(1:K, collapse=")(")
		partition <- paste("(", partition, ")", sep="")
		return(partition)
	}
	pairs <- matrix(nrow=0, ncol=2)
	for(i in 1:length(sep)){
		if(sep[i]==0) pairs <- rbind(pairs, idx_to_pair(i, K))
	}	

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
			j <- length(groups)
			if(!found) groups[[j+1]] <- p
		}
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
		return("Unk")
	}
	return(partition)
}
	

partition_to_sep <- function(partition, K){
	grps <-unlist(strsplit(partition, split="[()]"))
	grps <- grps[!grps==""]
	sep <- rep(1, K*(K-1)/2)
	for(j in 1:length(grps)){
		grp <- as.numeric(unlist(strsplit(grps[j], split=" ")))
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
