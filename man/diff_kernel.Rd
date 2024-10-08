% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff_kernel.R
\name{diff_kernel}
\alias{diff_kernel}
\title{diff_kernel}
\usage{
diff_kernel(W, nu, ld)
}
\arguments{
\item{W}{adjacency matrix (or Matrix)}

\item{nu}{strictly positive start measure on the vertices}

\item{ld}{2^ld is count of markov chain steps}
}
\value{
diffusion kernel as m(M)atrix
}
\description{
A function, given adjacency matrix, start measure, and observation time,
returns the corresponding diffusion kernel.
}
