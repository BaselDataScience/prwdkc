#' construct adjacency matrix from data
#'
#' @param dat data to cluster, either a dataframe or matrix (rows form the observations)
#' @param K count of nearest neighbors, defaults to floor(log(nrow(dat)))
#'
#' @return sparse adjacency matrix
#' @export
#'
#' @examples
#' W <- data2adj(seeds)
data2adj <- function(dat, K = floor(log(nrow(dat)))) {
  Matrix::sparseMatrix(i=rep(seq.int(nrow(dat)),K),
                       j=FNN::get.knn(dat, k=K)$nn.index,
                       x=1)
}
