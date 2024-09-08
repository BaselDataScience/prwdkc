#' find the optimal clustering
#'
#' @param X data matrix, rows are the vertices to be clustered
#' @param C matrix, each column is one clustering result
#'
#' @return optimal clustering vector according to the Calinski-Harabasz index
#'
#' @examples
optcluster <- function(X, C) {
  ch <- apply(C, 2, function(Y) fpc::calinhara(x=C, clustering = Y))
  C[,which.max(ch)]
}
