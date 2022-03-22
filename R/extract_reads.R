# A function to extract the reads containing a given k-mer
# and its reverse complement from bam files

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

# A function that generates a Manhattan plot of the results given k-mer positions, p-values, and a .fai index file
# sam_df: A sam-like data.frame of reads to which significant k-mers matched, as returned by the function add_pvalues
# fai_file: A .fai index file that is used to read chromosome sizes from
# pattern: a regular expression used to restrict the chromosomes that will be displayed on the Manhattan plot
kmer_manhattan <- function(sam_df, fai_file, pattern = NULL, signals = NULL) {
	fai_info <- parse_fai(fai_file, pattern)

	# Filtering sam_df be removing records that don't match pattern
	sam_df <- sam_df[grepl(pattern, sam_df$RNAME), ]

	# Getting the plotting position from the reference name and position
	sam_df$plot_pos <- fai_info$start[sam_df$RNAME] + sam_df$ref_pos
	sam_df$log10_p <- -log10(sam_df$pvalue)

	if(!is.null(signals) && nrow(signals)) {
		signals$first_plot <- fai_info$start[signals$first_chrom] + signals$first_pos
		signals$last_plot <- fai_info$start[signals$last_chrom] + signals$last_pos
		signals$max_plot <- fai_info$start[signals$max_chrom] + signals$max_pos
	}

	# Plotting the results
	output <- ggplot2::ggplot(sam_df) +
		ggplot2::geom_point(mapping = ggplot2::aes(x = plot_pos, y = log10_p, color = MAPQ)) +
		ggplot2::scale_x_continuous(name = "Chromosome",
					    breaks = fai_info$label_pos, 
					    labels = names(fai_info$start),
					    limits = c(0, sum(fai_info$lengths))) +
		ggplot2::ylab("-log10(p-value)") +
		ggplot2::theme_bw()

	if(!is.null(signals) && nrow(signals)) {
		output <- output + ggplot2::geom_vline(data = signals, aes(xintercept = max_plot))
		output <- output + ggplot2::geom_rect(data = signals,
						      ymin = 0, ymax = max(sam_df$log10_p),
						      aes(xmin = first_plot, xmax = last_plot),
						      color = "red", alpha = 0.4)
	}

	return(output)
}

# A function to read the output of GAPIT
read_gapit <- function(filename) {
	read.table(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE,
	colClasses = c("NULL", "integer", "numeric", "numeric", "NULL", "NULL", "NULL", "NULL", "numeric", "NULL"))
}

gapit_manhattan <- function(gapit_data, fai_file, pattern = NULL, signals = NULL, single_panel = TRUE, threshold = 5) {

	fai_info <- parse_fai(fai_file, pattern)

	# Getting the plotting position from the reference name and position
	gapit_data$first_chrom <- paste0("Gm", ifelse(gapit_data$Chromosome < 10, "0", ""), gapit_data$Chromosome)
	gapit_data$plot_pos <- fai_info$start[gapit_data$first_chrom] + gapit_data$Position
	gapit_data$log10_p <- -log10(gapit_data$FDR_Adjusted_P.values)

	if(!is.null(signals) && nrow(signals)) {
		signals$first_plot <- fai_info$start[signals$first_chrom] + signals$first_pos
		signals$last_plot <- fai_info$start[signals$last_chrom] + signals$last_pos
		signals$max_plot <- fai_info$start[signals$max_chrom] + signals$max_pos
	}

	# Setting a variable to alternate the colors on chromosomes
	gapit_data$colvar <- ifelse(gapit_data$Chromosome %% 2 == 0, "even", "odd")

	# If single_panel is TRUE, then all the data is plotted along a single axis
	if(single_panel) {
		output <- ggplot2::ggplot(gapit_data) +
			geom_point(mapping = aes(x = plot_pos, y = log10_p, color = colvar)) +
			scale_x_continuous(name = "Chromosome",
					   breaks = fai_info$label_pos,
					   labels = names(fai_info$start)) +
			ylab("-log10(p-value)") +
			geom_hline(yintercept = 5) +
			theme_bw() +
			theme(panel.grid.minor = element_blank())

		if(!is.null(signals) && nrow(signals)) {
			output <- output + ggplot2::geom_vline(data = signals, aes(xintercept = max_plot), linetype = 2)
		}

	} else {
	# If single_panel is not TRUE, there is one panel per chromosome
		output <- ggplot2::ggplot(gapit_data) +
			geom_point(mapping = aes(x = Position, y = log10_p)) +
			facet_wrap(~first_chrom, ncol = 5) +
			scale_x_continuous(name = "Position (bp)") +
			ylab("-log10(p-value)") +
			geom_hline(yintercept = threshold, linetype = 2) +
			theme_bw() +
			theme(panel.grid.minor = element_blank())

		if(!is.null(signals) && nrow(signals)) {
			output <- output + ggplot2::geom_vline(data = signals, aes(xintercept = max_plot))
			output <- output + ggplot2::geom_rect(data = signals,
							      ymin = 0, ymax = max(gapit_data$log10_p),
							      aes(xmin = first_pos, xmax = last_pos),
							      color = "red", alpha = 0.4)
		}
	}

	return(output)
}

# A function that takes a data.frame of GWAS results and returns a data.frame of signals that were found


