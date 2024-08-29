# find diffusion time for cluster separation
td_opt <- function(X, W, k, J) {
  P <- W / rowSums(W)  # Transition matrix
  xi <- colSums(P)  # ξ = νTP

  V <- matrix(nrow = J+1, ncol = nrow(W))
  P_1 <- diag(1/(xi+1)) %*% (P + t(P))
  for (j in 1:(J+1)) {
    V[j,] <- stats::kmeans(P_1, k=k)$cluster
    P_1 <- P_1 %*% P_1
  }
  
  ch <- apply(V, 1, function(Y) fpc::calinhara(x=X, clustering = Y))
  jmax <- which.max(ch)
  return(2^(jmax-1))
}
