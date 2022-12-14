% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kmer_plot.R
\name{format_haplotypes}
\alias{format_haplotypes}
\title{Format haplotypes for plotting using k-mer overlap information}
\usage{
format_haplotypes(haplotypes, overlaps)
}
\arguments{
\item{haplotypes}{A character vector of haplotypes to prepare for
plotting.}

\item{overlaps}{A list of IRanges object indicating the positions
of overlapping k-mers of interest and their associated -log10(p-value),
as returned by \code{\link{match_kmers}}. Note that a column called
"log10p", containing the -log10(p-value), is mandatory for format_haplotypes
to work, whereas it is optional for \code{\link{match_kmers}}. The names
of the elements of this list of IRanges must correspond to the haplotypes,
as the overlaps object is effectively queried by name.}
}
\value{
A list of data.frames with as many elements as there are haplotypes.
  Each such data.frame contains three columns:
  \itemize{
    \item pos: the position of the nucleotide along the haplotype
    \item nuc: the nucleotide located at that position
    \item log10p: the p-value associated with the most significant k-mer
                  overlapping that position. If no significant k-mer
                  overlaps that position, then the value is \code{NA}.
  }
}
\description{
This function prepares a set of haplotypes for plotting using the
function \code{\link{grid.haplotypes}} by formatting each haplotype
as a data.frame with one row per nucleotide and the associated p-value
of that nucleotide.
}
\examples{
NULL
}
