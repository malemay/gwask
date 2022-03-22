# A function that takes a data.frame returned by the add_pvalues function
# and returns a standard data.frame for Manhattan plotting
format_kmer_gwas <- function(sam_df) {

}

# A function that takes a fai index and returns a list of data that contains information about
# chromosome lengths so it can be used for Manhattan plotting
# fai_file: the fai index file to parse
# pattern: a regular expression to use to keep only chromosomes matching a given pattern
#          if NULL, all chromosomes are kept
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

