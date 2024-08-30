#' find optimal diffusion time for cluster separation
#'
#' @param X data matrix, rows representing the vertices to be clustered
#' @param W adjacency matrix
#' @param k number of clusters
#' @param J maximum power of 2 to consider
#'
#' @return optimal power of 2 to observe good cluster separation
#' @export
#'
#' @examples
optstop <- function(X, W, k, J) {
  X <- as.matrix(X)
  # Input checks
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix or dataframe")
  }

  # Input checks
  if ((!is.matrix(W) && !methods::is(W, 'Matrix')) || (!is.numeric(W) && !(is.numeric(W@x))) || any(W < 0)) {
    stop("W must be a numeric matrix with non-negative entries")
  }
  
  if (!is.numeric(k) || length(k) != 1 || k <= 0 || k != round(k)) {
    stop("k must be a positive integer")
  }
  
  if (!is.numeric(J) || length(J) != 1 || J <= 0 || J != round(J)) {
    stop("J must be a positive integer")
  }
  
  P <- W / rowSums(W)  # Transition matrix
  xi <- colSums(P)  # ξ = νTP

  # store clustering results for squares of P_1
  V <- matrix(nrow = J+1, ncol = nrow(W))
  P_1 <- diag(1/(xi+1)) %*% (P + t(P))
  for (j in 1:(J+1)) {
    V[j,] <- stats::kmeans(P_1, centers=k)$cluster
    P_1 <- P_1 %*% P_1
  }
  
  ch <- apply(V, 1, function(Y) fpc::calinhara(x=X, clustering = Y))
  jmax <- which.max(ch)
  return(jmax-1)
}
