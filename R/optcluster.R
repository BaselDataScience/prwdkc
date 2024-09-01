optcluster <- function(C, X) {
  ch <- apply(C, 2, function(Y) fpc::calinhara(x=C, clustering = Y))
  C[,which.max(ch)]
}
