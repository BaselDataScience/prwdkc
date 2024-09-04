#' prwdk clustering function
#'
#' @param W adjacency matrix (or Matrix)
#' @param k number of clusters to determine
#' @param nu start distribution on the vertices
#' @param ld 2^ld is stop time for determining clusters
#'
#' @return clustering, a vector indicating per vertex the assigned cluster number
#' @export
#'
#' @examples
prwdkc <- function(W, k, nu, ld) {
  N <- nrow(W)
  if (length(nu)==1) nu <- rep(nu, N)

  # Input checks
  if ((!is.matrix(W) && !methods::is(W, 'Matrix')) || (!is.numeric(W) && !(is.numeric(W@x))) || any(W < 0)) {
    stop("W must be a numeric matrix with non-negative entries")
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

  if (k > N) {
    stop("k cannot be larger than the number of vertices in the graph")
  }

  if (k > N / 2) {
    warning("k is more than half the number of vertices. This may lead to overfitting.")
  }

  # Main function body
  # Step 1: Compute the parametrized random walk operator P(ν)
  P <- W / Matrix::rowSums(W)  # Transition matrix
  xi <- as.vector(Matrix::crossprod(nu, P))
  P_nu <- diag(1/(1 + xi/nu)) %*% (P + Matrix::tcrossprod(diag(1/nu), P) %*% diag(nu))

  # Step 2: Compute the parametrized random walk diffusion kernel
  while (ld > 0) {
    P_nu <- P_nu %*% P_nu
    ld <- ld-1
  }
  K_td_nu <- P_nu %*% diag(1/(nu + xi))

  # Step 3: Apply k-means clustering to the rows of K(ld,ν)
  kmeans_result <- stats::kmeans(K_td_nu, centers = k, iter.max = 25)

  # Step 4: Return the k-partition from the clustering
  return(kmeans_result$cluster)
}
