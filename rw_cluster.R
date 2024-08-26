graph_clustering <- function(W, k, nu, td) {
  # Step 1: Compute the parametrized random walk operator P(ν)
  P <- W / rowSums(W)  # Transition matrix
  xi <- as.vector(t(nu) %*% P)  # ξ = νTP

  P_nu <- solve(diag(nrow(W)) + diag(xi) / nu) %*% (P + diag(1/nu) %*% t(P) %*% diag(nu))
  
  # Step 2: Compute the parametrized random walk diffusion kernel K(td,ν)
  K_td_nu <- expm::expm(td * log(P_nu)) %*% solve(diag(nu + xi))
  
  # Step 3: Apply k-means clustering to the rows of K(td,ν)
  kmeans_result <- kmeans(K_td_nu, centers = k)
  
  # Step 4: Obtain the k-partition Vtd based on the clustering result
  Vtd <- kmeans_result$cluster
  
  # Step 5: Return the partition
  return(Vtd)
}