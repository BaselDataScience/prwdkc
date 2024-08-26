# 8 node example

P <- matrix(c(rep(0,7), 0.6, 1, 0, 0.4, 1, rep(0,8), 1, rep(0,10), 0.4, rep(0,6), 0.8, rep(0,7), 0.2,
              rep(0,3), 0.6, 0, 0, 1, rep(0,3), 1, rep(0,6)), nrow=8)
det(P)

# symmetrize and make it Markov chain again
Psym <- (P+t(P))/2
Q <- diag(1/rowSums(Psym)) %*% Psym
x <- eigen(t(Q))$vectors[,1]
all.equal(y <- x/sum(x), as.vector(y %*% Q)) # y is stationary distribution

# stationary distribution p
p <- c(1, solve(t(P[-1,-1])-diag(nrow(P)-1), -P[1, -1]))
