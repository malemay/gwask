% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ld_functions.R
\name{cluster_ld}
\alias{cluster_ld}
\title{Greedy clustering of an LD matrix}
\usage{
cluster_ld(ld_matrix)
}
\arguments{
\item{ld_matrix}{A matrix of pairwise LD values between a set of k-mers}
}
\value{
A numeric matrix with the pairwise LD values computed between all
  k-mers, clustered such that k-mers with high LD are grouped together.
}
\description{
This function has been obsoleted by a much faster C implementation
}
\examples{
NULL
}
