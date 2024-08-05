graph_clustering <- function(W, k, nu, td) {
  N <- nrow(W)
  
  # Step 1: Compute the parametrized random walk operator P(ν)
  D <- diag(rowSums(W))
  P_nu <- solve(D) %*% W %*% diag(nu)
  
  # Step 2: Compute the parametrized random walk diffusion kernel K(td,ν)
  K_td_nu <- expm::expm(td * (P_nu - diag(N)))
  
  # Step 3: Apply k-means clustering to the rows of K(td,ν)
  kmeans_result <- kmeans(K_td_nu, centers = k)
  
  # Step 4: Obtain the k-partition Vtd based on the clustering result
  Vtd <- kmeans_result$cluster
  
  # Step 5: Return the partition
  return(Vtd)
}