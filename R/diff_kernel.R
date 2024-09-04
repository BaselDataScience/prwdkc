#' diff_kernel
#'
#' @description
#' A function, given adjacency matrix, start measure, and observation time,
#' returns the corresponding diffusion kernel.
#'
#' @param W adjacency matrix (or Matrix)
#' @param nu strictly positive start measure on the vertices
#' @param ld 2^ld is count of markov chain steps
#'
#' @return diffusion kernel as m(M)atrix
#' @export
#'
#' @examples
diff_kernel <- function(W, nu, ld) {
  # Input checks
  if ((!is.matrix(W) && !methods::is(W, 'Matrix')) || (!is.numeric(W) && !(is.numeric(W@x))) || any(W < 0)) {
    stop("W must be a numeric matrix with non-negative entries")
  }

  if (!is.numeric(nu) || !length(nu) %in% c(1,nrow(W)) || any(nu <= 0)) {
    stop("nu must be a positive numeric vector with length 1 or length equal to the number of rows in W")
  }

  if (!is.numeric(ld) || length(ld) != 1 || ld < 0 || ld != round(ld)) {
    stop("ld must be a nonnegative integer")
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
}
