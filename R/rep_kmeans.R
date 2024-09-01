#' rep_kmeans
#'
#' Repeats kmeans clusterings on a diffusion kernel matrix.
#'
#' @param Kd Diffusion kernel matrix.
#' @param k Number of clusters to determine.
#' @param n Number of repeated kmeans clusterings, default is 100
#'
#' @return numeric matrix of n cluster classifications.
#'   Each clustering is represented by one column of the matrix.
#' @export
#'
#' @examples
rep_kmeans <- function(Kd, k, n=100) {
  res <- vapply(seq.int(n), function(X) stats::kmeans(Kd, centers=k)$cluster, integer(nrow(Kd)))
}
