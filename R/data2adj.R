#' construct adjacency matrix from data
#'
#' @param dat data to cluster, either a dataframe or matrix (rows form the observations)
#'
#' @return sparse adjacency matrix
#' @export
#'
#' @examples
#' W <- data2adj(seeds)
data2adj <- function(dat) {
  N <- nrow(dat)
  k <- floor(log(nrow(dat)))
  Matrix::sparseMatrix(i=rep(seq.int(N),k),
                       j=FNN::get.knn(dat, k=k)$nn.index,
                       x=1)
}
