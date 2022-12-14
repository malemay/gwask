% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kmer_plot.R
\name{link_phenotypes}
\alias{link_phenotypes}
\title{Link haplotypes to their observed phenotypes}
\usage{
link_phenotypes(sequences, haplotypes, phenotypes, id_column, phenotype_column)
}
\arguments{
\item{sequences}{A named character vector of sequences observed
in a set of samples, as returned by \code{\link{read_fasta}}.
The names of the elements in the vector must correspond to
samples that are found in the phenotypes data.frame.}

\item{haplotypes}{A vector of haplotypes of interest observed
in the dataset. There must be a perfect correspondence between
these haplotypes and the sequences observed in the dataset.}

\item{phenotypes}{A data.frame of phenotypes observed in a set
of samples matching the ones for which sequences are supplied.
Must minimally include columns corresponding to id_column and
phenotype_column.}

\item{id_column}{A character of length one. The name of the column
in phenotypes that contains the sample IDs corresponding to the
names of the elements in \code{sequences}.}

\item{phenotype_column}{Similar to id_column, but for the column
containing the trait to use.}
}
\value{
A data.frame containing four columns:
\itemize{
   \item sample: the name of the sample
   \item haplotype: the sequence of the haplotype for that sample
   \item haplotype_id: a numeric identifier for the haplotype.
                       If NA, then the haplotype did not match
                       any of the reference haplotypes.
   \item phenotype: the value of the phenotype observed for the
                    trait of interest.
}
}
\description{
This function takes sequences observed in a set of samples, a set
of haplotypes of interest, and phenotypic data pertaining to the
samples, and returns a data.frame that links the haplotype
and its observed phenotype.
}
\examples{
NULL
}
