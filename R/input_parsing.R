#' Read the output of GAPIT
#'
#' This function reads the output of an analysis with GAPIT.
#' Only relevant fields (Chromosome, Position, P.value, FDR_Adjusted_P-values)
#' are kept to increase reading efficiency.
#'
#' @param filename Character of length one. The name of the .csv file
#'   containing the results of the GAPIT analysis.
#'
#' @return a data.frame containing the values read from the file
#'
#' @export
#' @examples
#' NULL
read_gapit <- function(filename) {
	read.table(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE,
	colClasses = c("NULL", "integer", "numeric", "numeric", "NULL", "NULL", "NULL", "NULL", "numeric", "NULL"))
}

#' Parse a .fai fasta index file
#'
#' This function takes a .fai index and returns some data that is useful
#' for generating Manhattan plots.
#'
#' @param fai_file Character of length one. The path to the .fai index file.
#' @param pattern Character of length one. A regular expression used to select
#'   the chromosomes or scaffolds to keep. If NULL, all entries in the .fai
#'   index are kept
#'
#' @return A list of length three with the following elements:
#'   \itemize{
#'     \item lengths: Numeric. The length of each sequence
#'     \item start: Numeric. The (cumulative) start position of that 
#'                  sequence along the genome
#'     \item label_pos: Numeric. The position to use for plotting chromosome
#'                      labels in the Manhattan plot.
#'   }
#'   All of the vectors in the list are named according to the chromosome/scaffold
#'
#' @export
#' @examples
#' NULL
parse_fai <- function(fai_file, pattern = NULL) {
	fai <- read.table(fai_file, header = FALSE, stringsAsFactors = FALSE,
			  col.names = c("NAME", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"))

	if(!is.null(pattern)) fai <- fai[grepl(pattern, fai[[1]]), ]
	
	# Extracting the starting positions of each chromosome from the cumulative sum of their lengths
	chrom_lengths <- fai$LENGTH

	chrom_start <- cumsum(chrom_lengths)
	chrom_start <- c(0, chrom_start[-(length(chrom_start))])

	label_positions <- chrom_start + 0.5 * chrom_lengths

	# All vectors have the names of the sequences to allow lookup of the info by name
	names(chrom_start) <- names(chrom_lengths)  <- names(label_positions) <- fai$NAME

	# Returning a list with useful information
	list(lengths = chrom_lengths, start = chrom_start, label_pos = label_positions)
}

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
	output <- makeGRangesFromDataFrame(sam_df,
					   start.field = "ref_pos",
					   seqnames.field = "RNAME",
					   end.field = "end",
					   keep.extra.columns = TRUE)

	# Setting the seqlevels and seqlengths of the output GRanges from the .fai index
	seqlevels(output) <- seqlevels(fai_info)
	seqlengths(output) <- seqlengths(fai_info)

	# Pruning the levels according to the pattern
	if(!is.null(pattern)) seqlevels(output) <- grep(pattern, seqlevels(output), value = TRUE)

	# Adding a column for the -log10(p-value)
	output$log10p <- -log10(output$pvalue)

	return(output)
}

#' Format GAPIT GWAS results for generating Manhattan plots
#'
#' This function takes a data.frame returned by the \code{\link{read_gapit}} function
#' and returns a standard data.frame for Manhattan plotting
#'
#' @param sam_df A data.frame of GAPIT analysis results as returned
#'   by \code{\link{read_gapit}}.
#' @param fai_file A character of length one. A .fai index file that
#'   will be queried with the \code{\link{parse_fai}} function to
#'   query information about the reference genome and translate
#'   positions on specific input sequences to positions along the
#'   whole reference.
#' @param chromosomes a character vector of the same length
#'   as the number of unique values in the Chromosome column
#'   of the GAPIT output. The integer chromosome values in the
#'   results returned by GAPIT will be used to index into that
#'   character vector to translate integer values into chromosome
#'   names.
#' @param pattern A character of length one. A regular expression
#'   used to match the name of the chromosome or scaffolds in the
#'   .fai index file. Matching sequences will be kept. If NULL
#'   (the default), all sequences are kept.
#'
#' @return A data.frame that contains the following standard columns
#'   and can be used for generating a Manhattan plot using the 
#'   function manhattan_plot.
#'   \itemize{
#'     \item manhattan_chrom Character. The name of the chromosome
#'     \item manhattan_rpos  Numeric. The position along the chromsome
#'       expressed as a cumulative position along the reference
#'     \item manhattan_cpos Numeric. The position on the reference
#'       expressed in terms of position along a particular chromosome
#'     \item mahattan_log10p The p-value used for plotting, obtained by
#'       computing -log10 of the p-value.
#'     \item manhattan_even Logical. \code{TRUE} if the chromosome should be
#'       considered even-valued, \code{FALSE} otherwise. Used to
#'       distinguish the chromosome that are next to each other by color.
#'   }
#'
#'   These data columns are guaranted to be the first ones, in that order,
#'   in the output data.frame. Any other columns that were present in the
#'   input are appended to the output data.frame as is.
#'
#' @export
#' @examples
#' NULL
format_gapit_gwas <- function(gapit_df, fai_file, chromosomes, pattern = NULL) {
	# Reading the .fai index
	fai_info <- parse_fai(fai_file, pattern)

	# Checking that the number of chromosomes in gapit df matches
	# the number of chromosomes in fai_info and the number of chromsomes
	# in the chromosomes character vector
	n_chrom <- length(unique(gapit_df$Chromosome))
	stopifnot(n_chrom == length(chromosomes) && n_chrom == length(fai_info$start))

	# Initializing an output data.frame to which the metadata columns will be appended
	output_df <- data.frame(manhattan_chrom = chromosomes[gapit_df$Chromosome],
				stringsAsFactors = FALSE)

	# Adding the other colums to the output data.frame
	output_df$manhattan_rpos <- fai_info$start[output_df$manhattan_chrom] +  gapit_df$Position
	output_df$manhattan_cpos <- gapit_df$Position
	output_df$manhattan_log10p <- -log10(gapit_df$FDR_Adjusted_P.values)
	output_df$manhattan_even   <- as.integer(factor(output_df$manhattan_chrom, levels = names(fai_info$start))) %% 2 == 0

	# Returning the output
	return(cbind(output_df, gapit_df))
}

