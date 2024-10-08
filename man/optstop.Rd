% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optstop.R
\name{optstop}
\alias{optstop}
\title{find optimal diffusion time for cluster separation}
\usage{
optstop(X, W, k, J)
}
\arguments{
\item{X}{data matrix, rows representing the vertices to be clustered}

\item{W}{adjacency matrix}

\item{k}{number of clusters}

\item{J}{maximum power of 2 to consider}
}
\value{
optimal power of 2 to observe good cluster separation
}
\description{
find optimal diffusion time for cluster separation
}
