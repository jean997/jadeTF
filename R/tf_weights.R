#' Solve trendfiltering problem with weights.
#' @description This function minimizes
#' \deqn{\frac{1}{2} \Vert w*(y - \theta) \Vert^{2}_{2} + \lambda\Vert D\theta\Vert_{1}}{
#' 1/2 || w(y - \theta) ||^2 + \lambda_||D\theta||_1}
#' Where D is a trendfiltering pendalty matrix.
#' This function borrows code heavily from the \code{genlasso} package and
#' also returns a \code{genlasso} object.
#'
#' @param y Vector of data.
#' @param pos Vector of positions.
#' @param wts Vector of weights
#' @param ord Order of polynomial to fit.
#' @param ... See \code{\link[genlasso]{trendfilter}} for other options.
#'
#' @return Returns a \code{genlasso} object with additional elements \code{weights}
#' with the weights provided.
#' @export
trendfilter_weights <- function(y, pos, wts=NULL, ord=0, approx=FALSE, maxsteps=2000,
                             minlam=0, rtol=1e-7, btol=1e-7, eps=1e-4, verbose=FALSE) {

  if (missing(y)) stop("y is missing.")
  if (!is.numeric(y)) stop("y must be numeric.")
  if (any(diff(pos)==0)) stop("pos must contain distinct values.")

  if(any(wts<=0)) stop("All wts must be positive.")
  if(is.null(wts)){
    #Equivalent to trendfilter with no X. This code is straight out of genlasso package
    n = length(y)
    if (is.null(pos)) D = genlasso::getDtfSparse(n,ord)
    else D = genlasso::getDtfPosSparse(n,ord,pos)
    out = genlasso:::dualpathWideSparse(y,D,NULL,approx,maxsteps,minlam,rtol,btol,verbose)

    # Compute beta, fit, y, bls (this would have been done by
    # the dualpath function)
    out$beta = as.matrix(y - t(D)%*%out$u)
    colnames(out$beta) = colnames(out$u)
    out$fit = out$beta
    out$y = y
    out$bls = y

    # Hijack the pathobjs component (to introduce what would
    # have been put here from the dualpath function)
    out$pathobjs$n0 = n
    out$pathobjs$y0 = y
    out$pathobjs$j = 1:n
    out$pathobjs$D0 = D
    out$pathobjs$coldif = 0
  }else{
    n = length(y)
    if (is.null(pos)) D = getDtfSparse(n,ord)
    else D = getDtfPosSparse(n,ord,pos)

    y2 = wts*y

    Xi = diag(1/wts)
    D2 = as(D %*% Xi, "sparseMatrix")

    out = genlasso:::dualpathWideSparse(y2,D2,NULL,approx,maxsteps,minlam,rtol,btol,verbose)

    # Compute beta, fit, y, bls (this would have been done by
    # the dualpath function)
    out$beta = Xi%*%as.matrix(y2 - t(D2)%*%out$u)
    colnames(out$beta) = colnames(out$u)
    out$fit = out$beta
    out$y = y
    out$bls = Xi %*% y2
    out$weights = wts

    # Hijack the pathobjs component (to introduce what would
    # have been put here from the dualpath function)
    out$pathobjs$n0 = n
    out$pathobjs$y0 = y2
    out$pathobjs$j = 1:n
    out$pathobjs$D0 = D2
    out$pathobjs$coldif = 0
  }
  out$pos=pos
  out$ord = ord
  class(out) = c("trendfilter", "genlasso", "list")
  return(out)
}
