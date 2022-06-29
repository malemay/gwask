#' Format k-mer GWAS results for downstream analyses
#'
#' This function takes a data.frame returned by the \code{\link{add_pvalues}}
#' function and returns a GRanges object suitable for downstream analyses 
#' and plotting.
#'
#' @param sam_df A data.frame of k-mer matches to alignments annotated
#'   with p-values, as returned by the function \code{\link{add_pvalues}}
#' @param ref_fasta A character of length one. A fasta file that
#'   will be queried with the \code{\link[Rsamtools]{scanFaIndex}} function to
#'   query information about the reference genome and set the seqinfo fields
#'   of the output GRanges object.
#' @param pattern A character of length one. A regular expression
#'   used to match seqlevels of the chromsomes or scaffolds in the
#'   fasta index. Matching sequences will be kept. If NULL
#'   (the default), all sequences are kept. In the current implementation,
#'   an error will be thrown if changing the seqlevels results in matches
#'   being pruned.
#' @param min_mapq A numeric of length one. The minimum mapping quality
#'   of the alignment on which the k-mer was found for it to be kept.
#'   The default value (0) implies no filtering.
#'
#' @return A GRanges object that contains the information on matching
#'   k-mers and can be used for downstream analyses including the
#'   function \code{\link{manhattan_plot}}. Must contain a field
#'   called log10_p that stores the -log10(p-value) of the k-mer.
#'
#' @export
#' @examples
#' NULL
format_kmer_gwas <- function(sam_df, ref_fasta, pattern = NULL, min_mapq = 0) {

	# Reading the .fai index info
	stopifnot(file.exists(paste0(ref_fasta, ".fai")))
	fai_info <- Rsamtools::scanFaIndex(ref_fasta)

	# Calculating the end of the ranges from the length of the kmers
	sam_df$end <- sam_df$ref_pos + nchar(sam_df$kmer) - 1

	# Filtering based on the MAPQ
	sam_df <- sam_df[sam_df$MAPQ >= min_mapq, ]

	# Coercing the input data.frame to a GRanges object
	output <- GenomicRanges::makeGRangesFromDataFrame(sam_df,
							  start.field = "ref_pos",
							  seqnames.field = "RNAME",
							  end.field = "end",
							  keep.extra.columns = TRUE)

	# Setting the seqlevels and seqlengths of the output GRanges from the .fai index
	GenomeInfoDb::seqlevels(output) <- GenomeInfoDb::seqlevels(fai_info)
	GenomeInfoDb::seqlengths(output) <- GenomeInfoDb::seqlengths(fai_info)

	# Pruning the levels according to the pattern
	if(!is.null(pattern)) GenomeInfoDb::seqlevels(output) <- grep(pattern, GenomeInfoDb::seqlevels(output), value = TRUE)

	# Adding a column for the -log10(p-value)
	output$log10p <- -log10(output$pvalue)

	return(output)
}

#' Format GAPIT GWAS results for generating Manhattan plots
#'
#' This function takes a CSV file of GAPIT GWAS results and formats
#' them as a GRanges object for processing with other functions, including
#' \code{\link{manhattan_plot}}.
#'
#' Only relevant fields (Chromosome, Position, P.value, FDR_Adjusted_P-values)
#' are kept to increase reading efficiency.
#'
#' @param filename Character of length one. The name of the .csv file
#'   containing the results of the GAPIT analysis.
#' @param ref_fasta A character of length one. A fasta file that
#'   will be queried with the \code{\link[Rsamtools]{scanFaIndex}} function to
#'   query information about the reference genome and set the seqinfo fields
#'   of the output GRanges object.
#' @param chromosomes a character vector of the same length
#'   as the number of unique values in the Chromosome column
#'   of the GAPIT output. The integer chromosome values in the
#'   results returned by GAPIT will be used to index into that
#'   character vector to translate integer values into chromosome
#'   names.
#' @param pattern A character of length one. A regular expression
#'   used to match seqlevels of the chromsomes or scaffolds in the
#'   fasta index. Matching sequences will be kept. If NULL
#'   (the default), all sequences are kept. In the current implementation,
#'   an error will be thrown if changing the seqlevels results in matches
#'   being pruned.
#'
#' @return A GRanges object that contains the information on associated
#'   markers and can be used for downstream analyses including the
#'   function \code{\link{manhattan_plot}}. Must contain a field
#'   called log10_p that stores the -log10(p-value) of the marker.
#'
#' @export
#' @examples
#' NULL
format_gapit_gwas <- function(filename, ref_fasta, chromosomes, pattern = NULL) {

	# Reading the file with the GAPIT results
	gapit_df <- read.table(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE,
			       colClasses = c("NULL", "integer", "numeric", "numeric", "NULL", "NULL", "NULL", "NULL", "numeric", "NULL"))

	# Reading the .fai index info
	stopifnot(file.exists(paste0(ref_fasta, ".fai")))
	fai_info <- Rsamtools::scanFaIndex(ref_fasta)

	# Checking that the number of chromosomes in gapit_df matches
	# the number of chromosomes in the chromosomes character vector
	n_chrom <- length(unique(gapit_df$Chromosome))
	stopifnot(n_chrom == length(chromosomes))

	# Translating the numeric values of the input file into chromosome names
	gapit_df$Chromosome <- chromosomes[gapit_df$Chromosome]

	# Coercing the gapit_df data.frame to a GRanges object
	output <- GenomicRanges::makeGRangesFromDataFrame(gapit_df,
							  start.field = "Position",
							  seqnames.field = "Chromosome",
							  end.field = "Position",
							  keep.extra.columns = TRUE)

	# Setting the seqlevels and seqlengths of the output GRanges from the .fai index
	GenomeInfoDb::seqlevels(output) <- GenomeInfoDb::seqlevels(fai_info)
	GenomeInfoDb::seqlengths(output) <- GenomeInfoDb::seqlengths(fai_info)

	# Pruning the levels according to the pattern
	if(!is.null(pattern)) GenomeInfoDb::seqlevels(output) <- grep(pattern, GenomeInfoDb::seqlevels(output), value = TRUE)

	# Adding a column for the -log10(corrected p-value) for downstream applications
	output$log10p <- -log10(gapit_df$FDR_Adjusted_P.values)

	# Returning the output
	return(output)
}

