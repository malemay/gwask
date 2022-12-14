% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ld_functions.R
\name{is_valid_ld}
\alias{is_valid_ld}
\title{Check the validity of an LD matrix}
\usage{
is_valid_ld(ld_matrix)
}
\arguments{
\item{ld_matrix}{A square numeric matrix with each cell indicating the pairwise
linkage disequilibrium between two k-mers.}
}
\value{
A logical value indicating whether the input LD matrix is valid.
}
\description{
This is a convenience function for verifying the validity of an LD matrix
used by any function that takes such a matrix as input. It is not exported
by the package and therefore only meant to be used internally.
}
\examples{
NULL
}
