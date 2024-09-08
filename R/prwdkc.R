#' prwdk clustering function
#'
#' @param X matrix or dataframe to be clustered, one element per row
#' @param k number of clusters to determine
#' @param nu start distribution on the vertices, defaulting to equidistribution
#' @param ld 2^ld is stop time for determining clusters, default computation by optstop()
#' @param K to compute K-NN matrix, default floor(log(nrow(X)))
#'
#' @return clustering, a vector indicating per vertex the assigned cluster number
#' @export
#'
#' @examples
prwdkc <- function(X, k, nu=1, ld=optstop(X, data2adj(X, K), k=k, J=15), K=floor(log(nrow(X)))) {
  N <- nrow(X)
  if (length(nu)==1) nu <- rep(nu, N)

  # Input checks
  if ((!is.data.frame(X) && !is.matrix(X) && !methods::is(X, 'Matrix')) ||
      (!is.numeric(X) && any(vapply(X, function(x) !is.numeric(x), logical(1))))) {
    stop("X must be a numeric matrix or dataframe")
  }

  if (!is.numeric(k) || length(k) != 1 || k <= 0 || k != round(k)) {
    stop("k must be a positive integer")
  }

  if (!is.numeric(nu) || length(nu) != N || any(nu <= 0)) {
    stop("nu must be a positive numeric vector with length 1 or length equal to the number of rows in W")
  }

  if (!is.numeric(ld) || length(ld) != 1 || ld < 0 || ld != round(ld)) {
    stop("ld must be a nonnegative integer")
  }

  if (!is.numeric(K) || length(K) != 1 || K <= 0 || K != round(K)) {
    stop("K must be a positive integer")
  }

  if (k > N) {
    stop("k cannot be larger than the number of vertices in the graph")
  }

  if (k > N / 2) {
    warning("k is more than half the number of vertices. This may lead to overfitting.")
  }

  # Main function body
  # Compute the adjacency matrix
  W <- data2adj(X, K)

  # Compute the parametrized random walk diffusion kernel
  K_td_nu <- diff_kernel(W=W, nu=nu, ld=ld)

  # Repeatedly apply k-means clustering to the rows of the kernel
  kmeans_result <- rep_kmeans(K_td_nu, k=k)

  # Return the best partition from the clustering
  return(optcluster(X, kmeans_result))
}
