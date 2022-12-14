% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{pvalue_tx_grob}
\alias{pvalue_tx_grob}
\title{Arrange a transcript grob and several p-value grobs in the same plot}
\usage{
pvalue_tx_grob(
  pvalue_grobs,
  xrange = NULL,
  xchrom = NULL,
  genes,
  transcripts,
  exons,
  cds,
  draw_arrows = FALSE,
  first_tx_only = FALSE,
  tx_fraction = 0.1,
  strand_colors = c("skyblue", "orange"),
  margins = c(5.1, 4.1, 4.1, 2.1)
)
}
\arguments{
\item{pvalue_grobs}{A list of p-value grobs that will be arrange in
viewport rows}

\item{xrange}{A GRanges object used to restrict the plotting region.
Typically, it will be a GWAS signal or a gene known to be involved
in the phenotype studied. If NULL (default), the scale is set automatically
to encompass all the values provided in the pvalue_grobs.}

\item{xchrom}{A character of length one indicating what chromosome is being
plotted. This parameter must be passed if xrange is NULL in order for
the GRanges object to use for subsetting genes to be set.}

\item{genes}{A GRanges object containing the known genes for the
reference genome being studied, as obtained from processing
a TxDb object of the GenomicFeatures package.}

\item{transcripts}{A GRangesList object containing the possible
transcripts for the set of genes in the reference genome being used.
Each element of the GRangesList object is a GRanges object containing
the known transcripts for a particular gene, with the name of the
list element being the name of the gene.}

\item{exons}{A GRangesList object containing the documented exons
in the reference genome being used. Each element of the list contains
the exons for a particular transcript, with the name of the list
element being the name of the transcript.}

\item{cds}{A GRangesList object containing the documented coding sequences
in the reference genome being used. Each element of the list contains
the coding sequences for a particular transcript, with the name of the list
element being the name of the transcript.}

\item{draw_arrows}{A logical indicating whether arrows should be drawn
in the middle of condig sequence regions to indicate the direction
of transcription.}

\item{first_tx_only}{A logical. Whether to plot only the first transcript
of a given gene (default: FALSE).}

\item{tx_fraction}{A numeric value between 0 and 1 representing the
proportion of the plot reserved for plotting the transcript(s).}

\item{strand_colors}{A vector of length 2 whose values can be interpreted
as colors. The first color will be used to color exons in genes that
are transcribed on the + strand, and the second for genes that are
transcribed on the - strand.}

\item{margins}{The margins used for the plotting viewport. Default value
is c(5.1, 4.1, 4.1, 2.1) (the default for \code{\link[grid]{plotViewport}})}
}
\value{
NULL, invisibly. The function is called for its side-effect of
  plotting.
}
\description{
Arrange a transcript grob and several p-value grobs in the same plot
}
\examples{
NULL
}
