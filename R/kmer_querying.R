# A function to extract the reads containing a given k-mer
# and its reverse complement from bam files

#' extract_reads
#'
#' To be completed
#'
#' @param kmer_sequence To be completed
#' @param bam_file To be completed
#' @param samtools_path To be completed
#' @param output To be completed
#'
#' @return To be completed
#'
#' @export
#' @examples
#' NULL
extract_reads <- function(kmer_sequence, bam_file, samtools_path, output) {
	# First extracting the reads with the kmer sequence in the forward orientation
	command <- paste0(samtools_path, " view ", bam_file, " | grep ", kmer_sequence, " > forward_", output)
	message("Running command: ", command)
	system(command)

	# Then running it on the reverse-complemented kmer sequence
	kmer_sequence <- revcomp(kmer_sequence)
	command <- paste0(samtools_path, " view ", bam_file, " | grep ", kmer_sequence, " > reverse_", output)
	message("Running command: ", command)
	system(command)
}


#' revcomp
#'
#' To be completed
#'
#' @param sequence To be completed
#'
#' @return To be completed
#'
#' @export
#' @examples
#' NULL
revcomp <- function(sequence) {

	# Checking that only one sequence is provided
	stopifnot(length(sequence) == 1)

	# The lookup table that will be used for replacement
	rep_table <- c("A" = "T",
		       "T" = "A",
		       "G" = "C",
		       "C" = "G",
		       "N" = "N")

	# Splitting the sequence into its constituent nucleotides
	sequence <- strsplit(sequence, "")[[1]]

	# Replacing the nucleotides by their complement
	sequence <- rep_table[sequence]

	# Returning the inverted sequence
	paste0(rev(sequence), collapse = "")
}

# A function that extracts the names of the samples containing a given k-mer
# pav_table: a file containing the presence/absence data of k-mers for the samples analyzed
# kmer_list: a file containing a list of k-mers (one per line) sorted by significant (from most to least significant)
# index: the index of the k-mer in kmer_list to extract
# 
# Return value: a list, with the first element being the kmer sequence and the second being a character vector of samples

#' get_samples
#'
#' To be completed
#'
#' @param pav_table To be completed
#' @param kmer_list To be completed
#' @param index To be completed
#'
#' @return To be completed
#'
#' @export
#' @examples
#' NULL
get_samples <- function(pav_table, kmer_list, index) {
	# First reading the PAV matrix and line names from pav_table
	# We need to read the sample names separately because they might get modified if we store them as column names
	pav_data <- read.table(pav_table, header = FALSE, skip = 1, stringsAsFactors = FALSE)
	line_names <- readLines(pav_table, n = 1)
	line_names <- strsplit(line_names, "\t")[[1]][-1]
	stopifnot(length(line_names) == (ncol(pav_data) - 1))

	# Then we read the set of k-mers
	kmers <- readLines(kmer_list)

	# We check that all the k-mers are in the PAV table, and vice-versa
	#stopifnot(all(kmers %in% pav_data[, 1]))
	#stopifnot(all(pav_data[, 1] %in% kmers))

	# Now we extract the k-mer that is of interest for this function call
	kmer  <- kmers[index]
	kmer_line <- which(pav_data[[1]] == kmer)
	stopifnot(length(kmer_line) == 1)

	# Extracting the samples that have the kmer of interest
	samples <- line_names[as.logical(pav_data[kmer_line, -1])]
	stopifnot(length(samples) > 0)

	return(list(kmer = kmer, samples = samples))
}

# A function that reads a pseudo-sam file (sam files without a header)
# and returns a data.frame containing the data
# Optional fields are stripped from the input
# filename: the psuedo-sam file to use as input

#' read_sam
#'
#' To be completed
#'
#' @param filename To be completed
#'
#' @return To be completed
#'
#' @export
#' @examples
#' NULL
read_sam <- function(filename) {

	# Handling the case where the file was empty
	if(!file.size(filename)) {
		output <- data.frame(QNAME = character(),
				     FLAG = integer(),
				     RNAME = character(),
				     POS = integer(),
				     MAPQ = integer(),
				     CIGAR = character(),
				     RNEXT = character(),
				     PNEXT = integer(),
				     TLEN = integer(),
				     SEQ = character(),
				     QUAL = character(),
				     stringsAsFactors = FALSE)
		return(output)
	}

	# We use col.names = paste0("V", 1:100) to force read.table to assume at least 100 input fields
	output <- read.table(filename, header = FALSE, 
			     stringsAsFactors = FALSE,
			     col.names = paste0("V", 1:100),
			     fill = TRUE)[, 1:11]

	colnames(output) <- c("QNAME", "FLAG", "RNAME", "POS",
			      "MAPQ", "CIGAR","RNEXT", "PNEXT",
			      "TLEN", "SEQ", "QUAL")

	output
}

# A function that takes the paths to a set of sam files and returns a data.frame with
# metadata about these files

#' make_sampath_df
#'
#' To be completed
#'
#' @param paths To be completed
#'
#' @return To be completed
#'
#' @export
#' @examples
#' NULL
make_sampath_df <- function(paths) {
	output = data.frame(path = paths,
			    kmer = sub("kmer_([ATGC]{31})_[0-9]+/.*", "\\1", paths),
			    sample = sub(".*(reverse|forward)_(.*)\\.sam", "\\2", paths),
			    orientation = sub(".*(reverse|forward).*", "\\1", paths),
			    stringsAsFactors = FALSE)
	output
}

# A function that takes a data.frame generated from the make_sampath_df function,
# reads the sam files, and returns a data.frame suitable for generating a Manhattan plot-like
# figure for the analysis of kmers with GWAS

#' parse_sam_data
#'
#' To be completed
#'
#' @param sam_paths To be completed
#'
#' @return To be completed
#'
#' @export
#' @examples
#' NULL
parse_sam_data <- function(sam_paths) {
	# Initializing a list that we will use to store the data.frames generated from each sam file
	sam_list <- list()

	# We initialize an empty data.frame for each sample/k-mer combination
	for(i in unique(paste0(sam_paths$sample, sam_paths$kmer))) {
		sam_list[[i]] <- data.frame(QNAME = character(),
					    FLAG = integer(),
					    RNAME = character(),
					    POS = integer(),
					    MAPQ = integer(),
					    CIGAR = character(),
					    RNEXT = character(),
					    PNEXT = integer(),
					    TLEN = integer(),
					    SEQ = character(),
					    QUAL = character(),
					    orientation = character(),
					    sample = character(),
					    kmer = character(),
					    stringsAsFactors = FALSE)
	}

	# Then we loop over all the input files to fill those data.frames
	for(i in 1:nrow(sam_paths)) {
		i_df <- read_sam(sam_paths[i, "path"])
		list_key <- paste0(sam_paths[i, "sample"], sam_paths[i, "kmer"]) 

		if(nrow(i_df)) {
			# The case if the sam file had at least one record
			i_df$orientation <- sam_paths[i, "orientation"]
			i_df$sample <- sam_paths[i, "sample"]
			i_df$kmer <- sam_paths[i, "kmer"]
		} else {
			# The case if the sam file was empty
			i_df$orientation <- character()
			i_df$sample <- character()
			i_df$kmer <- character()
		}

		# We add the data.frame at the appropriate location in sam_list
		sam_list[[list_key]] <- rbind(sam_list[[list_key]], i_df)
	}

	# All sample/kmer combinations should have at least one record
	stopifnot(all(sapply(sam_list, nrow) > 0))

	# We loop over the elements of sam_list to determine the major orientation for each kmer
	for(i in 1:length(sam_list)) {
		if(sum(sam_list[[i]]$orientation == "forward") >= sum(sam_list[[i]]$orientation == "reverse")) {
			sam_list[[i]]$major <- "forward"
		} else {
			sam_list[[i]]$major <- "reverse"
		}
	}

	# We genrate a data.frame by "rbinding" all the observations together
	sam_df <- do.call("rbind", sam_list)

	# We classify all the records according to whether they match the major or minor orientation
	sam_df$is_major <- ifelse(sam_df$orientation == sam_df$major, "major", "minor")

	# We also add a column that represents the k-mer sequence as it was found in the read (forward or reverse)
	sam_df$matching_seq <- character(nrow(sam_df))

	for(i in 1:nrow(sam_df)) {
		if(sam_df[i, "orientation"] == "forward") {
			sam_df[i, "matching_seq"] <- sam_df[i, "kmer"]
		} else {
			sam_df[i, "matching_seq"] <- revcomp(sam_df[i, "kmer"])
		}
	}

	return(sam_df)
}

# A function that takes the raw information returned from parse_sam_data and filters
# it be finding the start position of the k-mers on the reference and keeping only
# unique positions. By default, it keeps the read with the highest MAPQ among
# all the ones with the same mapping position

#' filter_sam_df
#'
#' To be completed
#'
#' @param sam_df To be completed
#'
#' @return To be completed
#'
#' @export
#' @examples
#' NULL
filter_sam_df <- function(sam_df) {
	# We use the sequence of the k-mer as it matched on the read to find its position in the reference
	# This should normalize the position for all reads containing the same k-mer
	sam_df$read_pos <- integer(nrow(sam_df))

	for(i in 1:nrow(sam_df)) {

		read_pos <- regexpr(sam_df[i, "matching_seq"], sam_df[i, "SEQ"], fixed = TRUE)
		read_pos <- as.numeric(read_pos)
		stopifnot(all(read_pos > 0))
		stopifnot(length(read_pos) == 1)

		sam_df[i, "read_pos"] <- read_pos
	}

	# The read_pos column is the position of the k-mer within the read
	# We use this information to get the position on the reference
	sam_df$ref_pos <- sam_df$POS + sam_df$read_pos - 1

	# We want to keep only kmers with a unique position in the dataset
	# For each such kmer we will keep the read with the largest mapping quality
	sam_df <- sam_df[order(sam_df$MAPQ, decreasing = TRUE), ]
	sam_df <- sam_df[!duplicated(sam_df[, c("RNAME", "ref_pos", "kmer")]), ]
	sam_df
}

# A function that adds p-values from the pvalue-sorted file of GWAS results to the input data.frame

#' add_pvalues
#'
#' To be completed
#'
#' @param sam_df To be completed
#' @param input To be completed
#'
#' @return To be completed
#'
#' @export
#' @examples
#' NULL
add_pvalues <- function(sam_df, input) {
	# Reading the input file
	pvalues <- read.table(input, header = FALSE, stringsAsFactors = FALSE)

	# Obtaining the k-mer sequence
	pvalues$kmer_seq <- sub("([ATGC]{31})_.*", "\\1", pvalues[[2]])

	# All k-mer sequences should be 31 characters
	stopifnot(all(nchar(pvalues$kmer_seq) == 31))

	# Obtaining the p-values for sam_df bny matching the k-mer sequences
	sam_df$pvalue <- pvalues[match(sam_df$kmer, pvalues$kmer_seq), 9]

	# We return a data.frame sorted by p-value for most significant to least significant
	sam_df[order(sam_df$pvalue), ]
}

