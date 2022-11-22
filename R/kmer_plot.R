#' Read a fasta file into a named character vector
#'
#' This function reads a fasta file into a named character vector.
#' Its functionality is limited as it is required that all sequences
#' be on a single line as the function assumes one line per sequence.
#'
#' @param input_fasta A character. The name of the fasta file to
#'   read the input from.
#'
#' @return A character vector with as many elements as there are
#'   sequences. The name of each element is the name of the sequence
#'   (the sequence of characters that follows ">").
#'
#' @export
#' @examples
#' NULL
read_fasta <- function(input_fasta) {
	fasta_lines <- readLines(input_fasta)

	# Extracting the names of the samples
	sample_names <- grep("^>", fasta_lines, value = TRUE)
	sample_names <- sub("^>", "", sample_names)

	# Extracting the sequences themselves from the fasta file
	sequences <- grep("^>", fasta_lines, value = TRUE, invert = TRUE)

	# Checking that the number of samples matches the number of sequences
	if(length(sample_names) != length(sequences)) {
		stop("The number of samples in ", input_fasta, " does not match the number of\n",
		     "sequences. read_fasta() does not support sequences that run over several lines.")
	}

	# The names of the samples are the names of the sequence vector
	names(sequences) <- sample_names

	return(sequences)
}

#' Read the (sorted) k-mer p-values from the output of a k-mer GWAS analysis
#'
#' This function reads the p-values associated with a set of k-mers output by
#' a k-mers GWAS analysis. The k-mers must have been sorted by increasing p-value
#' beforehand as it supports reading only the first top "n" p-values instead
#' of the whole file. This function does check that the values read in are sorted
#' by increasing p-value, but this only guarantees that the first n p-values were
#' sorted.
#'
#' @param kmer_file A character. The name of the text file containing the k-mers
#'   and p-values. The k-mer IDs (which will be processed to extract only the k-mer
#'   sequence) must be int he second column, whereas the p-values must be in ninth
#'   (and last) column.
#' @param max_kmers The maximum number of k-mers to read from the file. The function
#'   will issue a warning if the data.frame returned contains fewer rows than the total
#'   number of records in the file.
#' @param kmer_length A numeric of length one. The length of the k-mers used in
#'   the analysis. Only used for sanity checks. The function will throw an error
#'   if any of the k-mers do not match this length.
#'
#' @return A data.frame with four columns:
#'   \itemize{
#'     \item{kmer: A character. The sequence of the k-mer.}
#'     \item{kmer_reverse: A character. The reverse sequence of the k-mer.}
#'     \item{pvalue: A numeric. The p-value associated with the k-mer.}
#'     \item{log10p: A numeric. the -log10 of the p-value associated with the k-mer.}
#'   }
#'
#' @export
#' @examples
#' NULL
read_kmer_pvalues <- function(kmer_file, max_kmers, kmer_length) {

	# Reading the k-mer sequences and p-values from file, using only relevant columns
	output <- read.table(kmer_file,
			     colClasses = c("NULL", "character", rep("NULL", 6), "numeric"),
			     nrows = max_kmers)

	# Checking if we have reached the maximum number of k-mers
	if(nrow(output) == max_kmers) {
		warning("read_kmer_ptavlues() reached maximum number of k-mers (", max_kmers, ") in ", kmer_file,
			"\nSome matching k-mers may be missed.")
	}

	# Naming the two columns
	names(output) <- c("kmer_id", "pvalue")

	# As we assume that the p-value column is sorted, we need to check for this
	if(is.unsorted(output$pvalue)) stop("ERROR: read_kmer_pvalues expects sorted p-values in ", kmer_file)

	# Extracting the k-mer sequence from the k-mer ID
	output$kmer <- sub("[^ATGC]+", "", output$kmer_id)

	# Getting the reverse complement of the k-mer sequence
	output$kmer_reverse <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(output$kmer)))

	# Checking that all k-mer sequences have the expected length
	if(!all(nchar(output$kmer) == kmer_length) || !all(nchar(output$kmer_reverse) == kmer_length)) {
		stop("Not all k-mers match expected length of ", kmer_length)
	}
	
	# Computing the -log10 of the p-value
	output$log10p <- -log10(output$pvalue)

	# Only keeping the necessary columns for returning the data.frame
	return(output[, c("kmer", "kmer_reverse", "pvalue", "log10p")])
}

#' Extract and filter the haplotypes from a set of sequences
#'
#' This function identifies the haplotypes (unique sequences) found
#' in a set of sequences (as read using \code{\link{read_fasta}})
#' and filters them based on their frequency.
#'
#' @param sequences A character vector of sequences. This vector is
#'   typically named, but need not be. A check is made that the
#'   sequences contain no other characters than those in the set ATGCN.
#' @param min_frequency A numeric of length one. The minimum number of
#'   times that a haplotype must occur for it to be kept. By default,
#'   no filterinf is done (min_frequency = 0).
#'
#' @return A character vector of haplotypes observed in the dataset with
#'   a frequency >= min_frequency.
#'
#' @export
#' @examples
#' NULL
get_haplotypes <- function(sequences, min_frequency) {

	# Check that the sequence only contains valid nucleotides (ATGCN)
	if(any(grepl("[^ATGCNatgcn]", sequences))) {
		stop("Error in get_haplotypes(): one or several sequences contain invalid character(s).")
	}

	# First compute the number of times that a given sequence occurs in the dataset
	haplotypes <- table(sequences)

	# Keeping only haplotypes that occur at least min_frequency times
	haplotypes <- haplotypes[haplotypes >= min_frequency]

	# Sorting the haplotypes by decreasing order of frequency before the clustering them
	haplotypes <- sort(haplotypes, decreasing = TRUE)
	haplotypes <- names(haplotypes) # We no longer care about the frequency, only about the sequence

	if(length(haplotypes) > 2) {
		haplotypes <- cluster_haplotypes(haplotypes)	
	}

	# We return the (possibly clustered) sequences
	return(haplotypes)
}

#' Cluster the most similar haplotypes together
#'
#' This function takes a character vector of strings and uses a
#' greedy clustering algorithm based on distance between strings
#' to cluster them together.
#'
#' @param sequences A character vector of sequences to cluster.
#'
#' @return A character vector with the same values as the input
#'   vector, but reordered to reflect similarity between strings.
#'
#' @examples
#' NULL
cluster_haplotypes <- function(sequences) {
	to_cluster <- sequences[-1]
	clustered <- sequences[1]

	while(length(to_cluster) > 0) {
		sdist <- adist(clustered[length(clustered)], to_cluster)
		c_index <- which.min(sdist)
		clustered <- c(clustered, to_cluster[c_index])
		to_cluster <- to_cluster[-c_index]
	}

	clustered
}

#' Link haplotypes to their observed phenotypes
#'
#' This function takes sequences observed in a set of samples, a set
#' of haplotypes of interest, and phenotypic data pertaining to the
#' samples, and returns a data.frame that links the haplotype
#' and its observed phenotype.
#'
#' @param sequences A named character vector of sequences observed
#'   in a set of samples, as returned by \code{\link{read_fasta}}.
#'   The names of the elements in the vector must correspond to
#'   samples that are found in the phenotypes data.frame.
#' @param haplotypes A vector of haplotypes of interest observed
#'   in the dataset. There must be a perfect correspondence between
#'   these haplotypes and the sequences observed in the dataset.
#' @param phenotypes A data.frame of phenotypes observed in a set
#'   of samples matching the ones for which sequences are supplied.
#'   Must minimally include columns corresponding to id_column and
#'   phenotype_column.
#' @param id_column A character of length one. The name of the column
#'   in phenotypes that contains the sample IDs corresponding to the
#'   names of the elements in \code{sequences}.
#' @param phenotype_column Similar to id_column, but for the column
#'   containing the trait to use.
#'
#' @return A data.frame containing four columns:
#' \itemize{
#'    \item sample: the name of the sample
#'    \item haplotype: the sequence of the haplotype for that sample
#'    \item haplotype_id: a numeric identifier for the haplotype.
#'                        If NA, then the haplotype did not match
#'                        any of the reference haplotypes.
#'    \item phenotype: the value of the phenotype observed for the
#'                     trait of interest.
#' }
#'
#' @export
#' @examples
#' NULL
link_phenotypes <- function(sequences, haplotypes, phenotypes, id_column, phenotype_column) {

	# Sanity checks
	if(!id_column %in% colnames(phenotypes)) stop("Column ", id_column, " not found in phenotypes.")
	if(!phenotype_column %in% colnames(phenotypes)) stop("Column ", phenotypes_column, " not found in phenotypes.")

	# Initialize the output data with the sample names and their haplotype
	haplotype_data <- data.frame(sample = names(sequences),
				     haplotype = unname(sequences),
				     stringsAsFactors = FALSE)

	# Assigning a haplotype to each sample
	haplotype_data$haplotype_id <- NA

	for(i in 1:nrow(haplotype_data)) {
		if(haplotype_data[i, "haplotype"] %in% haplotypes) {
			haplotype_data[i, "haplotype_id"] <- which(haplotypes == haplotype_data[i, "haplotype"])
		}
	}

	# Adding the phenotype for each sample
	haplotype_data$phenotype <- phenotypes[match(haplotype_data$sample, phenotypes[[id_column]]), phenotype_column]

	# Returning the data.frame
	return(haplotype_data)
}

#' Match a set of k-mers to positions on sequences
#'
#' This function find the positions of the matches of a set
#' of k-mers in a vector of sequences, optionally carrying
#' over metadata columns present in the k-mers dataset.
#'
#' @param sequences A character vector of DNA sequences to search
#'   for matches among the set of k-mers of interest.
#' @param kmers A data.frame of k-mers of interest, as returned by
#'   \code{\link{read_kmer_pvalues}}. Must minimally contain a column
#'   called "kmer" for the k-mer sequence and another called "kmer_reverse"
#'   for the reverse-complemented sequence of the k-mer, plus any column
#'   listed in \code{data_columns}.
#' @param kmer_length A numeric of length one. The length of the k-mers
#'   being considered. Used for sanity checks and for setting the width
#'   of the output IRanges.
#' @param data_columns A character vector of data colums to carry over
#'   from the \code{kmers} data.frame to the output IRanges object.
#'   If \code{NULL} (the default), then no data columns are carried over.
#'
#' @return A list of IRanges objects with as many elements as the number
#'   of sequences being queried. The name of the list elements corresponds
#'   to the sequence queried. Each IRanges object contains the positions
#'   where matching k-mers were identified. By default, no other information
#'   than this position is returned; data columns from the \code{kmers} object
#'   can optionally be carried over to the output through the \code{data_columns}
#'   argument.
#'
#' @export
#' @examples
#' NULL
match_kmers <- function(sequences, kmers, kmer_length, data_columns = NULL) {

	# Check for the appropriate format of the kmers dataset
	if(!all(c("kmer", "kmer_reverse") %in% colnames(kmers))) {
		stop("Error in match_kmers(): 'kmers' lacks the appropriate columns 'kmer' and/or 'kmer_reverse'")
	}

	# Check that kmer_length does match the length of the k-mers
	if(!all(nchar(kmers$kmer) == kmer_length) || !all(kmers$reverse_kmer)) {
		stop("Error in match_kmers(): the length of the kmers does not match the expected value of ", kmer_length)
	}

	# Creating a list that will have as many elements as there are sequences
	# This list will be output by the function
	overlap_list <- list()

	for(i in sequences) {
		overlap_list[[i]] <- kmers
		overlap_list[[i]]$fmatch   <- sapply(kmers$kmer, function(kmer) regexpr(kmer, i, fixed = TRUE))
		overlap_list[[i]]$rmatch   <- sapply(kmers$kmer_reverse, function(kmer) regexpr(kmer, i, fixed = TRUE)) 
		overlap_list[[i]]$matchpos <- pmax(overlap_list[[i]]$fmatch, overlap_list[[i]]$rmatch)
	}

	# Keeping only the k-mers for which there an overlap and relevant columns
	columns <- "matchpos"

	if(!is.null(data_columns)) {
		stopifnot(is.character(data_columns) && all(nchar(data_columns) > 0) && all(data_columns %in% colnames(kmers)))
		columns <- c(columns, data_columns)
	}

	overlap_list <- lapply(overlap_list, function(x) x[x$matchpos != -1, columns])

	# Coercing to an IRanges object and optionally appending data columns
	overlap_list <- lapply(overlap_list, function(x, data_columns) {
				       output <- IRanges::IRanges(start = x$matchpos, width = kmer_length)

				       # Optionally adding metaata columns (if any)
				       if(!is.null(data_columns)) {
					       S4Vectors::mcols(output) <- x[, data_columns]
					       names(S4Vectors::mcols(output)) <- data_columns
				       }

				       return(output)},
			     data_columns = data_columns)

	return(overlap_list)
}

#' Format haplotypes for plotting using k-mer overlap information
#'
#' This function prepares a set of haplotypes for plotting using the
#' function \code{\link{grid.haplotypes}} by formatting each haplotype
#' as a data.frame with one row per nucleotide and the associated p-value
#' of that nucleotide.
#'
#' @param haplotypes A character vector of haplotypes to prepare for
#'   plotting.
#' @param overlaps A list of IRanges object indicating the positions
#'   of overlapping k-mers of interest and their associated -log10(p-value),
#'   as returned by \code{\link{match_kmers}}. Note that a column called
#'   "log10p", containing the -log10(p-value), is mandatory for format_haplotypes
#'   to work, whereas it is optional for \code{\link{match_kmers}}. The names
#'   of the elements of this list of IRanges must correspond to the haplotypes,
#'   as the overlaps object is effectively queried by name.
#'
#' @return A list of data.frames with as many elements as there are haplotypes.
#'   Each such data.frame contains three columns:
#'   \itemize{
#'     \item pos: the position of the nucleotide along the haplotype
#'     \item nuc: the nucleotide located at that position
#'     \item log10p: the p-value associated with the most significant k-mer
#'                   overlapping that position. If no significant k-mer
#'                   overlaps that position, then the value is \code{NA}.
#'   }
#' @export
#' @examples
#' NULL
format_haplotypes <- function(haplotypes, overlaps) {

	# Testing if the log10p column is found in the overlaps IRanges objects
	if(!all(sapply(overlaps, function(x) "log10p" %in% names(S4Vectors::mcols(x))))) {
		stop("Error in format_haplotypes(): 'log10p' must be among the metadata columns of 'overlaps'")
	}

	# Testing if all haplotypes are represented in the overlaps dataset
	if(!all(haplotypes %in% names(overlaps))) {
		stop("Error in format_haplotypes(): not all haplotypes are represented in 'overlaps'")
	}

	# Iterate over all haplotypes using lapply
	output <- lapply(haplotypes, function(x, kmer_overlaps) {

				 # Initializing the output data.frame with each nucleotide and its position
				 output_df <- data.frame(pos = 1:nchar(x),
							 nuc = strsplit(x, "")[[1]],
							 log10p = NA,
							 stringsAsFactors = FALSE)

				 # Getting the maximum -log10(p-value) for that nucleotide (if any)
				 for(i in 1:nrow(output_df)) {
					 i_range <- IRanges::IRanges(start = output_df[i, "pos"], width = 1)
					 overlapped_pos <- IRanges::subsetByOverlaps(kmer_overlaps[[x]], i_range)
					 if(length(overlapped_pos)) {
						 output_df[i, "log10p"] <- max(S4Vectors::mcols(overlapped_pos)$log10p)
					 }
				 }

				 output_df
			     },
			     kmer_overlaps = overlaps)

	# Returning the list of data frames thus generated
	return(output)
}

#' Align a set of sequences using mafft
#'
#' This function takes a set of nucleotide sequences and aligns them
#' using the mafft program for multiple sequence alignment. The mafft
#' program is launched through the \code{system} interface as it is not
#' an \code{R} program, however the results are imported directly into
#' \code{R} and formatted for output.
#'
#' @param fasta_path A character. The path to the fasta file to use as
#'   an input file to mafft. This file is only temporarily used, but it
#'   is not deleted after the function executes (might be changed later).
#' @param sequences A character vector of nucleotide sequences to align.
#' @param mafft_path A character. The path to the \code{mafft} executable.
#' @param mafft_options A character. A set of options to pass to \code{mafft}
#'   between the name of the executable and the name of the fasta file containing
#'   the sequences to align. By default, this argument is an empty string.
#'
#' @return A character vector of aligned sequences which will all have the
#'   same length due to gaps introduced by the alignment procedure.
#'
#' @export
#' @examples
#' NULL
mafft_align <- function(fasta_path, sequences, mafft_path, mafft_options = "") {

	# Opening the fasta file that will be used as input for mafft
	fasta <- file(fasta_path, open = "w+")
	on.exit(close(fasta))

	# Writing the sequences in fasta format to fasta_path
	for(i in 1:length(sequences)) {
		cat(paste0(">seq", i, "\n"), file = fasta)
		cat(sequences[i], "\n", file = fasta)
	}

	# Prepare the mafft command and launch it using the system interface
	# The output of mafft is read directly into an R variable
	mafft_command <- paste0(mafft_path, " ", mafft_options, " ", fasta_path)
	alignment <- system(mafft_command, intern = TRUE)

	# Sanity check that mafft has not re-ordered the sequences
	alignment_num <- regmatches(alignment, regexpr("[0-9]+", alignment))
	stopifnot(all(diff(as.numeric(alignment_num)) == 1))

	# Formatting the alignment for output
	alignment <- strsplit(paste0(alignment, collapse = ""), split = ">seq[0-9]+")[[1]]
	alignment <- alignment[nchar(alignment) > 0]

	# Checking that all alignments have indeed the same length
	stopifnot(length(unique(nchar(alignment))) == 1)

	return(alignment)
}

#' Adjust plotting position according to alignment gaps
#'
#' This function accepts a dataset of haplotypes formatted for plotting
#' (with \code{\link{grid.haplotypes}}), as returned by the function
#' \code{\link{format_haplotypes}}, as well as a multiple sequence
#' alignment produced for those haplotypes using \code{\link{mafft_align}}.
#' This information is used to update the formatted haplotypes by introducing
#' gaps (shown as dashes, '-') in the alignment prior to plotting.
#'
#' Once deleted positions have been identified, gaps are filled using the
#' function \code{\link{fill_gaps}}, which is not exported by the package.
#'
#' @param hapdata A list of data.frames containing information about haplotypes
#'   to plot, as returned by \code{\link{format_haplotypes}}. Each such data.frame
#'   must include the columns 'pos', 'nuc' and 'log10p'.
#' @param alignment A character vector of aliigned sequences to be used for adjusting
#'   the plotting positions and adding gaps in hapdata. The sequences must have a
#'   one-to-one correspondence with the ones in hapdata; currently no sanity check
#'   is implemented to make sure that this is indeed the case.
#'
#' @return A list of data.frames formatted as the input hapdata object, but with
#'   positions modified to take the alignment into account and gaps added.
#'   Gaps will be shown as a '-' nucleotide and \code{NA} as a -log10(p-value).
#'
#' @export
#' @examples
#' NULL
adjust_gaps <- function(hapdata, alignment) {

	# Some sanity checks on the hapdata object
	if(!all(sapply(hapdata, function(x) all(c("nuc", "pos", "log10p") %in% colnames(x))))) {
		stop("Error in adjust_gaps(): one or several required data columns missing in hapdata.")
	}

	# Identifying the positions of the gaps from the alignment and formatting to an IRanges object
	gaps <- stringr::str_locate_all(alignment, "-")
	gaps <- lapply(gaps, function(x) IRanges::IRanges(start = x[, 1], end = x[, 2]))
	gaps <- lapply(gaps, function(x) IRanges::reduce(x))

	# There should be as many gaps as there are haplotypes
	stopifnot(length(hapdata) == length(gaps))

	# We iterate over the haplotypes/gaps
	for(i in 1:length(hapdata)) {

		# We check if there are any gaps for this haplotype
		if(length(i_gaps <- gaps[[i]])) {

			# We iterate over the gaps that were found for this haplotype
			for(j in length(i_gaps)) {
				indices <- which(hapdata[[i]]$pos >= IRanges::start(i_gaps[j]))
				hapdata[[i]][indices, "pos"] <- hapdata[[i]][indices, "pos"] + IRanges::width(i_gaps[j])
			}
		}
	}

	# Now that the positions have been adjusted, we can fill the deleted positions with dashes
	hapdata <- fill_gaps(hapdata)

	return(hapdata)
}

#' Fill the deleted positions with dashes
#'
#' This function takes a list of data.frames that have been
#' adjusted for position by taking alignment gaps into account,
#' and fills those positions with gaps by setting the nucleotide
#' to a dash ('-') and the -log10(p-value) to \code{NA}. This
#' function is not exported for direct use, but is instead called
#' internally by \code{\link{adjust_gaps}} to simplify its code.
#'
#' @param hapdata A list of data.frames containing data about
#'   haplotypes to be plotted, possibly with missing positions
#'   due to alignment gaps.
#'
#' @return A list of data.frames similar in format to hapdata,
#'   in which positions that correspond to gaps have been filled.
#'
#' @examples
#' NULL
fill_gaps <- function(hapdata) {
	# We need to get the set of positions covered from 1 to the maximum across all haplotypes
	maxpos <- max(do.call("rbind", hapdata)$pos)

	# Then we loop over all the data.frames and fill the missing positions with dashes
	for(i in 1:length(hapdata)) {
		indices <- which(! 1:maxpos %in% hapdata[[i]]$pos)

		if(!length(indices)) next

		new_rows <- data.frame(pos = indices,
				       nuc = "-",
				       log10p = NA,
				       stringsAsFactors = FALSE)

		hapdata[[i]] <- rbind(hapdata[[i]], new_rows)
		hapdata[[i]] <- hapdata[[i]][order(hapdata[[i]]$pos), ]
	}

	return(hapdata)
}

#' Find positions that differ between two sequences
#'
#' This function takes as input a list of haplotypes that have been
#' prepared for plotting and identifies the positions. The positions
#' of the nucleotides must have been prepared with \code{\link{adjust_gaps}}
#' such that there is a one-to-one correspondence between positions.
#' Gaps are not flagged as differing between sequences.
#'
#' @param hapdata A list of data.frames prepared for plotting,
#'   as returned by \code{\link{adjust_gaps}}.
#'
#' @return A list of numeric vectors indicating the positions in
#'   hapdata where top consecutive haplotypes differ. The length
#'   of the returned list is therefore length(hapdata) - 1.
#'
#' @export
#' @examples
#' NULL
nucdiff <- function(hapdata) {

	# Checking that there are at least two haplotypes to compare
	if(length(hapdata) < 2) stop("nucdiff needs at least two haplotypes to compare")

	# Checking that all haplotypes have the same number of positions
	if(!length(unique(sapply(hapdata, nrow))) == 1) {
		stop("Error in nucdiff(): all input haplotypes must have the same length.")
	}

	# Creating a list that will contain the positions that differ
	output <- list()

	for(i in 1:(length(hapdata) - 1)) {
		if(!all(hapdata[[i]]$pos == hapdata[[i + 1]]$pos)) stop("All positions must be shared between haplotypes")

		# Identifying positions where a haplotype and the next one differ (without being a gap)
		output[[i]] <- which(hapdata[[i]]$nuc != hapdata[[i + 1]]$nuc & 
				     hapdata[[i]]$nuc != "-" &
				     hapdata[[i + 1]]$nuc != "-")
	}

	return(output)
}

#' Plot a set of haplotypes using grid functions
#'
#' This function uses grid functions to plot haplotype sequences at a
#' given locus. Each sequence is plotted in its own row viewport and
#' positions where consecutive haplotypes differ are optionally marked
#' through the \code{difflist} argument. Nucleotides are color-coded
#' according to the -log10(p-value) of the most significant k-mer
#' that they are part of, thus highlighting the positions of highly
#' significant k-mers and suggesting causal variants underlying
#' the presence/absence of those k-mers. Color mapping is controlled
#' through arguments that are passed to \code{\link{map_color}},
#' a function that is not exported by this package.
#'
#' @param hapdata A list of data.frames containing information about haplotypes
#'   to plot, as returned by \code{\link{adjust_gaps}}. All haplotypes must
#'   have the same number of positions.
#' @param difflist A list of positions where sequences differ between two
#'   consecutive haplotypes, as returned by \code{\link{nucdiff}}. The length
#'   of this list must be length(hapdata) - 1. This argument is used to
#'   represent sequence difference graphically. If \code{NULL} (the default),
#'   then sequence differences are not represented graphically.
#' @param fontsize The font size
#' @param position A character of length one indicating the chromosome and
#'   nucleotide position x scale, used for plotting the x-axis. This parameter
#'   must be specified as in samtools, using the format Chr:start-end. If NULL,
#'   then no x-axis is plotted.
#' @inheritParams map_color
#'
#' @return \code{NULL}, invisibly. This function is invoked for its plotting
#'   side-effects.
#'
#' @export
#' @examples
#' NULL
grid.haplotypes <- function(hapdata, difflist = NULL, fontsize = 8,
			    n_colors = 5, pal = "YlOrRd", scale_extend = 2,
			    position = NULL) {

	# Checking that all haplotypes have the same number of positions
	if(!length(unique(sapply(hapdata, nrow))) == 1) {
		stop("Error in grid.haplotypes(): all input haplotypes must have the same length.")
	}

	# Getting the number of haplotypes to plot
	n_hap <- length(hapdata)

	# Reformatting the data as a single data.frame
	for(i in 1:n_hap) hapdata[[i]]$hapnum <- i
	hapdata <- do.call("rbind", hapdata)

	# Mapping the -log10(pvalues) onto a color scale
	map_color_output <- map_color(values = hapdata$log10p,
				      pal = pal,
				      n_colors = n_colors + 2,
				      scale_extend = scale_extend,
				      n_colors_remove = 2)

	hapdata$color <- map_color_output$mapped_colors

	# Extracting informations for the x-axis (if applicable)
	if(!is.null(position)) {
		stopifnot(length(position) == 1 && is.character(position))
		chrom <- sub(":.*", "", position)
		stopifnot(nchar(chrom) > 0)

		xscale <- sub(".*:", "", position)
		xscale <- as.numeric(strsplit(xscale, "-")[[1]])
	}

	# Dividing the viewport in 9 sub-viewports in a 3x3 layout parts:
	# - The upper part will be used to plot the color scale
	# - The middle part acts as a buffer between the upper and lower viewports
	# - The lower part will be used to plot the haplotype sequences
	# - Three columns are also used to leave some space off to the left and right of the plotting region
	grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 3,
								     ncol = 3,
								     widths = grid::unit(c(5, 0.86, 1), c("null", "npc", "null")),
								     heights = grid::unit(c(0.5, 0.4, 0.1), "null"))))

	# Further dividing the top viewport into rows and columns
	# ----- Haplotype rows are interleaved with rows used to show sequence differences between consecutive haplotypes
	# ----- Hence the number of rows being 2 * n_hap - 1
	# ----- The number of columns is the maximum position across all haplotypes
	grid::pushViewport(grid::viewport(layout.pos.row = 1,
					  layout.pos.col = 2,
					  layout = grid::grid.layout(nrow = n_hap * 2 - 1,
								     ncol = max(hapdata$pos),
								     heights = grid::unit(rep(c(2, 1), length.out = n_hap * 2 - 1), "null")),
					  xscale = if(is.null(position)) c(0, 1) else xscale))

	# Adding the x-axis if applicable
	if(!is.null(position)) {
		grid::grid.xaxis()
		grid.text(paste0("Position along ", chrom, " (bp)"),
			  y = grid::unit(-3, "lines"))
	}

	# Iterating over the rows in hapdata to plot the nucleotides in the appropriate row and column
	for(i in 1:nrow(hapdata)) {
		grid::grid.text(hapdata[i, "nuc"],
				gp = grid::gpar(fontsize = fontsize,
						fontfamily = "mono",
						fontface = "bold",
						col = hapdata[i, "color"]),
				vp = grid::viewport(layout.pos.row = hapdata[i, "hapnum"] * 2 - 1,
						    layout.pos.col = hapdata[i, "pos"]))
	}

	# Adding haplotype labels in the left margin
	for(i in 1:n_hap) {
		grid::grid.text(paste0("Hap ", i),
				x = grid::unit(-0.02, "npc"),
				hjust = 1,
				vp = grid::viewport(layout.pos.row = i * 2 - 1))

	}

	# Also plotting the positions where the nucleotides differ (if difflist is provided)
	if(!is.null(difflist)) {

		# Check that the length of difflist is n_hap -1
		if(n_hap - length(difflist) != 1) {
			stop("Error in grid.haplotypes(): the length of difflist is not (number of haplotypes) - 1")
		}

		for(i in 1:length(difflist)) {
			if(!length(difflist[[i]])) next

			# Looping over the differing positions between haplotype i and the next one
			for(j in difflist[[i]]) {
				grid::grid.lines(x = 0.5,
						 y = c(0, 1),
						 vp = grid::viewport(layout.pos.row = i * 2, layout.pos.col = j))
			}
		}
	}

	# Finally plotting the color scale in the top viewport
	grid::upViewport() # Back in the viewport that is divided in two rows
	grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 2))
	grid.colorscale(breaks = map_color_output$breaks,
			base_palette = map_color_output$base_palette,
			label_text = expression(-log[10](italic(p))),
			round_digits = 0,
			direction = "horizontal")
	grid::upViewport()

	# Moving back into the top viewport
	grid::upViewport()

	return(invisible(NULL))
}

#' Plot a contingency table of observed phenotypes and haplotypes
#'
#' This function uses grid function to plot the co-occurence of phenotypes
#' and haplotypes. This can be used to analyze whether specific haplotypes
#' are linked to the occurence of certain phenotypes.
#'
#' @param phenodata A data.frame of observed haplotypes and the corresponding
#'   phenotype observed for a given trait in a set of sample. The data.frame must
#'   minimally contain a "phenoytpe" column with the phenotypic value and a
#'   "haplotype_id" column with a numeric identifier for the haplotype. Each
#'   row represents different sample.
#'
#' @return \code{NULL}, invisibly. This function is invoked for its plotting
#'   side-effects.
#'
#' @export
#' @examples
#' NULL
grid.phenotable <- function(phenodata) {

	# Checking for required columns in phenodata
	if(!all(c("phenotype", "haplotype_id") %in% colnames(phenodata))) {
		stop("Error in grid.phenotable(): missing columns in phenodata.")
	}

	# Computing a table of the phenotypes observed per haplotype
	phenotable <- table(phenodata$phenotype, phenodata$haplotype_id)

	# Creating a viewport with the required cells to print the data
	grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = nrow(phenotable) + 1,
								     ncol = ncol(phenotable) + 1,
								     widths =  grid::unit(c(0.2, rep(1, ncol(phenotable))), c("npc", rep("null", ncol(phenotable)))),
								     heights = grid::unit(c(0.2, rep(1, nrow(phenotable))), c("npc", rep("null", nrow(phenotable)))))))

	# Printing the IDs of the phenotypes
	for(i in 1:ncol(phenotable)) {
		grid::grid.text(paste0("Hap ", colnames(phenotable)[i]), vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i + 1))
	}

	# Printing the phenotype values themselves
	for(i in 1:nrow(phenotable)) {
		grid::grid.text(rownames(phenotable)[i],
				x = grid::unit(0.1, "npc"),
				hjust = 0,
				vp = grid::viewport(layout.pos.row = i + 1, layout.pos.col = 1))
	}

	# Printing the counts
	for(i in 1:nrow(phenotable)) {
		for(j in 1:ncol(phenotable)) {
			grid::grid.text(as.character(phenotable[i, j]),
					vp = grid::viewport(layout.pos.row = i + 1, layout.pos.col = j + 1))
		}
	}

	# Adding some lines to make the contingency table more interesting to look at
	for(i in 1:ncol(phenotable)) {
		grid::grid.lines(x = grid::unit(c(0, 0), "npc"),
				 y = grid::unit(c(0, 1), "npc"),
				 vp = grid::viewport(layout.pos.col = i + 1))
	}

	for(i in 1:nrow(phenotable)) {
		grid::grid.lines(x = grid::unit(c(0, 1), "npc"),
				 y = grid::unit(c(1, 1), "npc"),
				 vp = grid::viewport(layout.pos.row = i + 1))
	}

	# Returning to the top viewport
	grid::upViewport()

	return(invisible(NULL))
}

#' Map a set of numeric values onto a color palette
#'
#' This function takes a set of numeric values and maps them on a color palette
#' provided by the RColorBrewer package. \code{NA} values are mapped to black
#' as a default color for missing values.
#'
#' @param values A numeric vector of values to be mapped on the color scale.
#'   \code{NA} values are accepted and mapped to the "black" color.
#' @param n_colors A numeric. The number of colors to use for the mapping p-values
#'   onto a color scale. The maximum value allowed may differ depending on the
#'   particular palette used. This represents the number of colors in the base
#'   palette prior to applying the n_colors_remove option.
#' @param pal A character. The palette to use for to -log10(p-values) color scale.
#'   must correspond to a palette in the RColorBrewer package.
#' @param scale_extend A numeric. The number of units (of -log10(p-value)) to extend
#'   the limits of the color scale by at both ends to make sure that extreme values
#'   are mapped onto the scale.
#' @param n_colors_remove A numeric. The number of colors to remove from the beginning
#'   of the color scale, to remove paler hues. If this option is used, the length of the
#'   final palette will be n_colors - n_colors_remove.
#'
#' @return A list of length three comprising the following elements:
#'   \itemize{
#'     \item mapped_colors: a character vector of colors mapped from the values input vector
#'     \item breaks: a numeric vector of break positions used for mapping
#'     \item base_palette: the colors that were used for mapping the numeric values
#'   }
#'
#' @examples
#' NULL
map_color <- function(values, pal, n_colors, scale_extend, n_colors_remove = 0) {
	base_palette <- RColorBrewer::brewer.pal(n = n_colors, name = pal)
	if(n_colors_remove > 0) {
		base_palette <- base_palette[-(1:n_colors_remove)]
		n_colors <- length(base_palette)
	}

	# Getting the breaks based on the minimum/maximum values and the scale_extend parameter
	breaks <- seq(min(values, na.rm = TRUE) - scale_extend,
		      max(values, na.rm = TRUE) + scale_extend,
		      length.out = n_colors + 1) 

	# A vector of colors to use for plotting; NA values are changed to black
	mapped_colors <- base_palette[as.integer(cut(values, breaks = breaks, include.lowest = TRUE))]
	mapped_colors[is.na(mapped_colors)] <- "black"

	# We return a list with the mapped values, color levels and base palette colors
	list(mapped_colors = mapped_colors,
	     breaks = breaks,
	     base_palette = base_palette)
}

#' Plot the color scale used in a haplotype plot
#'
#' This function takes a set of numeric breaks and a set of colors
#' and uses this information to produce a graphical scale linking
#' those breaks to the colors used in a plot. This function is not
#' exported from the package and is therefore meant to be used only
#' internally by the package.
#'
#' @param breaks A numeric vector of breaks used for the color scale.
#'   There must be one more break than the number of colors
#' @param base_palette A character vector of colors used in mapping
#'   the numeric values onto the color scale.
#' @param label_text A character to use as a value for the axis label.
#' @param digits A numeric value indicating the number of digits
#'   that the breaks in the scale should be rounded to.
#' @param direction A character indicating whether the direction in
#'   which the color scale should be plotted. Supported values are
#'   "horizontal" and "vertical".
#' @param fontsize A numeric. The font size to use for the axis label
#'   and tick labels.
#'
#' @return \code{NULL}, invisibly. This function is invoked for its
#'   plotting side-effects.
#'
#' @examples
#' NULL
grid.colorscale <- function(breaks, base_palette, label_text, round_digits = 2,
			    direction = "horizontal", fontsize = 12) {

	# Checking that the number of breaks is one more than the number of colors
	if(length(breaks) != length(base_palette) + 1) {
		stop("Error in grid.colorscale(): length(breaks) should be equal to length(base_palette) + 1")
	}

	# Checking the value of direction
	if(!direction %in% c("horizontal", "vertical")) {
		stop("Error in grid.colorscale: direction must be one of 'horizontal' or 'vertical'")
	}

	horiz <- direction == "horizontal"

	# Getting the x-scale from the minimum and maximum values of the break labels
	xscale <- range(breaks)

	# Pushing a viewport with the appropriate scale
	grid::pushViewport(grid::viewport(width = 0.5,
					  height = 0.5,
					  xscale = xscale,
					  yscale = xscale))

	# Draw the scale colors as rectangles
	grid::grid.rect(x = if(horiz) breaks[-length(breaks)] else grid::unit(0, "npc"),
			width = if(horiz) diff(breaks) else grid::unit(1, "npc"),
			y = if(horiz) grid::unit(0, "npc") else breaks[-length(breaks)],
			height = if(horiz) grid::unit(1, "npc") else diff(breaks),
			default.units = "native",
			just = c(0, 0),
			gp = grid::gpar(fill = base_palette))

	# Adding the breaks for the scale
	if(horiz) {
		grid::grid.xaxis(at = breaks,
				 label = as.character(round(breaks, digits = round_digits)),
				 gp = grid::gpar(fontsize = fontsize))
	} else {
		grid::grid.yaxis(at = breaks,
				 label = as.character(round(breaks, digits = round_digits)),
				 gp = grid::gpar(fontsize = fontsize),
				 main = FALSE)
	}

	# Adding a name for the scale
	if(horiz) {
		grid::grid.text(label_text,
				y = grid::unit(-3, "lines"),
				gp = grid::gpar(fontsize = fontsize))
	} else {
		grid::grid.text(label_text,
				x = grid::unit(5, "lines"),
				gp = grid::gpar(fontsize = fontsize),
				rot = 90)
	}

	grid::upViewport()

	return(invisible(NULL))
}

