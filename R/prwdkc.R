prwdkc <- function(W, k, nu, td) {
  # Input checks
  if (!is.matrix(W) || !is.numeric(W) || any(W < 0)) {
    stop("W must be a numeric matrix with non-negative entries")
  }
  
  if (!is.integer(k) || length(k) != 1 || k <= 0) {
    stop("k must be a positive integer")
  }
  
  if (!is.numeric(nu) || length(nu) != nrow(W) || any(nu <= 0)) {
    stop("nu must be a positive numeric vector with length equal to the number of rows in W")
  }
  
  if (!is.integer(td) || length(td) != 1 || td < 0) {
    stop("td must be a nonnegative integer")
  }
  
  if (k > nrow(W)) {
    stop("k cannot be larger than the number of vertices in the graph")
  }
  
  if (k > nrow(W) / 2) {
    warning("k is more than half the number of vertices. This may lead to overfitting.")
  }
  
  # Main function body
  # Step 1: Compute the parametrized random walk operator P(ν)
  P <- W / rowSums(W)  # Transition matrix
  xi <- as.vector(t(nu) %*% P)  # ξ = νTP
  P_nu <- diag(1/(1 + xi/nu)) %*% (P + diag(1/nu) %*% t(P) %*% diag(nu))
  
  # Step 2: Compute the parametrized random walk diffusion kernel K(td,ν)
  K_td_nu <- (P_nu %^% td) %*% diag(1/(nu + xi))
  
  # Step 3: Apply k-means clustering to the rows of K(td,ν)
  kmeans_result <- kmeans(K_td_nu, centers = k)
  
  # Step 4: Obtain the k-partition Vtd based on the clustering result
  Vtd <- kmeans_result$cluster
  
  # Step 5: Return the partition
  return(Vtd)
}