# read seeds data

library(FNN)
library(igraph)
library(Matrix)
library(readr)

seeds <- readr::read_table('data/seeds_dataset.txt',
                           col_names = c('area',
                                         'perimeter',
                                         'compactness',
                                         'length',
                                         'width',
                                         'asymmetry',
                                         'groove',
                                         'variety'
                                         )
                           )

# construct KNN adjacency matrix W
N <- nrow(seeds)
K <- floor(log(N))
W <- Matrix::sparseMatrix(rep(1:N, K),
                          FNN::get.knn(seeds[,-8], k=K)$nn.index,
                          x=1)   # can be improved by using similarities instead of 1

