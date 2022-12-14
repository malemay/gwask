% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{manhattanGrob}
\alias{manhattanGrob}
\title{A function that returns a graphical object (grob) representing a manhattan plot}
\usage{
manhattanGrob(
  gwas_results,
  threshold = NULL,
  min_log10p = 0,
  ref_signals = NULL,
  new_signals = NULL,
  label_offset = 10^7,
  label_hjust = 0,
  numeric_chrom = FALSE,
  point_colors = c("blue", "red"),
  yexpand = c(0.03, 0.08),
  xexpand = c(0.03, 0.03),
  cex.points = 0.5,
  cex.lab = 1,
  margins = c(5.1, 4.1, 4.1, 2.1)
)
}
\arguments{
\item{gwas_results}{A GRanges object formatted with either \code{\link{format_kmer_gwas}}
or \code{\link{format_gapit_gwas}}. These functions
return GRanges containing standard columns and seqinfo in addition to GWAS type-specific
columns, therefore allowing a common interface for all GWAS types.}

\item{threshold}{Numeric. A single value indicating where the -log10(p-value)
threshold line should be plotted.}

\item{min_log10p}{Numeric. The minimum -log10(p-value) to display on the plot
can be used to filter out some non-significant markers for plotting speedup.}

\item{ref_signals}{A GenomicRanges object containing reference signals to be added
to the plot. If NULL (the default), then no signals are plotted.}

\item{new_signals}{A GenomicRanges object containing signals that were discovered
for this trait. When provided, this colors the labels of the reference signals
that were found in green, and those that were not found in red. Supplying
new_signals without supplying ref_signals throws an error.}

\item{label_offset}{A numeric value indicating the offset of the signal labels
in base pairs relative to the location of the signal itself. Used to add
signal labels without overlapping the dotted line itself.}

\item{label_hjust}{A numeric value indicating the hjust parameter for the
signal labels. Should be set according to the value used for label_offset.}

\item{numeric_chrom}{Logical. Whether chromosome names should be stripped from
their alphabetical component and converted to numeric values (default: FALSE).}

\item{point_colors}{A vector of length two whose values can be interpreted as colors.
the first value is used for the color of odd-numbered chromsomes, and the second
value for even-numbered chromosomes.}

\item{yexpand}{A numeric of length two. A fractional value relative to the
range of the y-axis used to expand the y-scale on either side. The first
value is the expansion factor at the bottom of the scale, and the second
value is the expansion factor at the top of the scale.}

\item{xexpand}{Same as yexpand but for the x-scale.}

\item{cex.points}{A numeric value. The expansion factor for the points.}

\item{cex.lab}{A numeric value. The expansion factor for the axis labels.}

\item{margins}{The margins used for the plotting viewport. Default value
is c(5.1, 4.1, 4.1, 2.1) (the default for \code{\link[grid]{plotViewport}})}
}
\value{
A gTree describing a Manhattan plot with the specified parameters.
}
\description{
A function that returns a graphical object (grob) representing a manhattan plot
}
\examples{
NULL
}
