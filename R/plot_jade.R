
plot_data <- function(counts, reads, pos, ylim=c(0, 1),
                      range=NULL, xlab="Position", ylab="Methylation Proportion"){
  if(!is.matrix(counts)){
				K=1
				p <- length(counts)
				counts <- matrix(counts, nrow=p)
				reads <- matrix(reads, nrow=p)
  }else{
		p <- dim(counts)[1]
		K <- dim(counts)[2]
  }
	cat(p, K, "\n")
  if(K==1) sv <- TRUE
  else sv <- FALSE

	cols=c("black", "red3", "blue2")
	if(K >3) cols <- c(cols, 4:(K-3))

	if(is.null(range)){
		range <- range(pos)
		idx <- 1:p
	}else{
		idx <- which(pos <= max(range) & pos >= min(range))
		pos <- pos[idx]
		cat(dim(fits), "\n")
	}
	p <- length(pos)
  data.full <- (counts/reads)[idx,]

	df <- data.frame(pos, data.full, reads)
  names(df) <- c("pos", paste("y", 1:K, sep=""), paste("r", 1:K, sep=""))
	h <- ggplot(df)
  for(j in 1:K){
			h <- h + geom_point(aes_string(x="pos", y=paste("y", j, sep=""),
			                               size=paste("r", j, sep="")), col=cols[j], alpha=1)
  }
  return(h + scale_size_area() + ylim(ylim)+ labs(y=ylab, x=xlab) + theme(legend.position="none",
                                                   text=element_text(size=20))
         )
}

#' Plot a JADE fit
#' @param fits A p by K matrix of fitted values (the \code{fits}) element of a JADE fit object.
#' @param pos Vector of positions.
#' @param y A p by K matrix of data. Optional. If provided, data points will be added to the plot.
#' Data can be plotted by itself using the \code{plot_data} function.
#' @param reads A p by K matrix of read counts or an indication of the relative amount of data for
#' each element of \code{y}. Line widths will be proportional to \code{reads}.
#' @param wsize Window size used to determine line widths.
#' @param sep.tab A table giving regions which are separated.
#' This can be obtained using \code{\link{get_separated_regions}}.
#' @return A \code{ggplot} object for the desred plot.
#' @export
plot_jade <- function(fits,  pos, y=NULL, reads=NULL, cols=NULL,
                       range=NULL, wsize=10, maxcov=Inf, take.log.cov=FALSE,
                       maxwidth=3, minwidth=0.5,  ylim=NULL, sep.tab=NULL,
                      ylab=NULL, xlab="Position", bg.color="blue"){

  if(!is.null(y)) plot.data=TRUE
	  else plot.data=FALSE

	if(class(fits)=="numeric"){
			p <- length(fits)
			K <- 1
			cat("Single variable.\n")
			if(!is.null(y)) y <- matrix(y, nrow=p)
			if(!is.null(reads)) reads <- matrix(reads, nrow=p)
			alpha <- 0.8
			sv<- TRUE
	}else{
		p <- dim(fits)[1]
		K <- dim(fits)[2]
		sv <- FALSE
		alpha <- 0.6
	}

	if(is.null(cols)){
	  cols=c("black", "red3", "blue2")
	  if(K >3) cols <- c(cols, 4:(K-3))
	}


	if(is.null(range)){
		range <- range(pos)
		idx <- 1:p
	}else{
		idx <- which(pos <= max(range) & pos >= min(range))
		fits <- fits[idx,, drop=FALSE]
		pos <- pos[idx]
		cat(dim(fits), "\n")
	}
	p <- length(pos)

	if(!is.null(y)) data.full <- y[idx,, drop=FALSE]
	   else data.full <- matrix(NA, nrow=p, ncol=K)

  if(!is.null(reads)) reads.full <- reads[idx,, drop=FALSE]
	  else reads.full <- matrix(1, nrow=p, ncol=K)

	if(is.null(ylim) & !plot.data) ylim=range(fits)
	  else if(is.null(ylim)) ylim=range(cbind(data.full, fits))

	window.coverage <- matrix(NA, nrow=p, ncol=K)
	for(i in 1:p){
		t <- pos[i]
		if(sum(pos >= t-(wsize/2) & pos <= t+(wsize/2)) ==1){
			 C <- reads.full[i, ]
		}else{
			C <- apply(reads.full[ pos >= t-(wsize/2) & pos <= t+(wsize/2),,drop=FALSE] , MARGIN=2, FUN=sum)
		}
		window.coverage[i, ] <- C
	}
	window.coverage[window.coverage > maxcov] <- maxcov
	window.coverage <- window.coverage + 2
	if(take.log.cov){
		window.coverage=log10(window.coverage)
		window.coverage[window.coverage== -Inf] <- 0
	}


	df <- data.frame(pos, fits, window.coverage, data.full, reads.full)
  names(df) <- c("pos", paste("t", 1:K, sep=""), paste("w", 1:K, sep=""), paste("y", 1:K, sep=""), paste("r", 1:K, sep=""))
	h <- ggplot(df)

	#Background colors
	if(!is.null(sep.tab)){
		h <- h + gg_sep_rect(sep.tab, color = bg.color, alpha=0.2)
	}

	for(j in 1:K){
		h <- h +  geom_line(aes_string(x="pos", y=paste("t", j, sep=""), size=paste("w", j, sep="")),
		                    col=cols[j], alpha=alpha)
	}

	if(plot.data){
	  for(j in 1:K){
	    h <- h + geom_point(aes_string(x="pos", y=paste("y", j, sep="")),
	                        col=cols[j], size=2, shape=1)
	  }
	}
	R <- h + scale_x_continuous(minor_breaks=pos) +
	  ylim(ylim) + scale_size(range=c(minwidth, maxwidth))
	if(is.null(xlab) & is.null(ylab)){
	  R <- R + theme(legend.position="none",
	                 axis.title.x=element_blank(),
	                 axis.title.y=element_blank(),
	                 text=element_text(size=20))
	}else if(is.null(xlab)){
	  R <- R + theme(legend.position="none",
	                 axis.title.x=element_blank(),
	                 text=element_text(size=20)) + labs(y=ylab)
	}else if(is.null(ylab)){
	  R <- R + theme(legend.position="none",
	                 axis.title.y=element_blank(),
	                 text=element_text(size=20)) + labs(x=xlab)
	}else{
	  R <- R + theme(legend.position="none",
	                 text=element_text(size=20)) + labs(x=xlab, y=ylab)
	}
	return(R)
}


gg_sep_rect <- function(sep.tab, color="blue", alpha=0.2){
  m.lower <- sep.tab$Start
  m.upper <- sep.tab$Stop
  M <- data.frame(m.lower, m.upper)
  names(M) <- c("m.lower", "m.upper")

 return(geom_rect(data=M, aes(xmin=m.lower,xmax=m.upper,ymin=-Inf,ymax=Inf),
                     fill=color, alpha=alpha))
}

