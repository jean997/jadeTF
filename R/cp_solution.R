#df <- clusterpath.l1.id(X)

find_solution_clusterpath_cpp <- function(df, lambda){
		K <- nlevels(df$row)
		p <- nlevels(df$col)
		beta.out <- rep(0, p*K)
		mybeta <- .C("read_cp_R",
		             alpha = as.double(df$alpha),
		             lam = as.double(df$lambda),
		             row = as.integer(df$row),
		             col = as.character(df$col),
		             gamma=as.double(lambda),
		             length = as.integer(nrow(df)),
		             k=as.integer(K),
		             beta=as.double(beta.out), package="jadeTF")
		beta <- matrix(mybeta$beta, nrow=p, byrow=TRUE)
		return(beta)
}
