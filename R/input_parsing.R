# A function that takes a data.frame returned by the add_pvalues function
# and returns a standard data.frame for Manhattan plotting
format_kmer_gwas <- function(sam_df) {

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

# A function to read the output of GAPIT
read_gapit <- function(filename) {
	read.table(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE,
	colClasses = c("NULL", "integer", "numeric", "numeric", "NULL", "NULL", "NULL", "NULL", "numeric", "NULL"))
}

# A function that takes a data.frame of GWAS results and returns a data.frame of signals that were found

