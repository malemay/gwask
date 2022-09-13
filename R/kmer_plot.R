#' Plot a table of phenotypes versus haplotypes
#'
#' Details
#'
#' @param haplotype_data To complete
#'
#' @return To complete
#'
#' @export
#' @examples
#' NULL
grid.phenotable <- function(haplotype_data) {
	# Computing a table of the phenotypes observed per haplotype
	phenotable <- table(haplotype_data$phenotype, haplotype_data$haplotype_id)

	# Creating a viewport with the required cells to print the data
	grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = nrow(phenotable) + 1,
								     ncol = ncol(phenotable) + 1)))

	# Printing the IDs of the phenotypes
	for(i in 1:ncol(phenotable)) {
		grid::grid.text(colnames(phenotable)[i], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i + 1))
	}

	# Printing the phenotype values themselves
	for(i in 1:nrow(phenotable)) {
		grid::grid.text(rownames(phenotable)[i], vp = grid::viewport(layout.pos.row = i + 1, layout.pos.col = 1))
	}

	# Printing the counts
	for(i in 1:nrow(phenotable)) {
		for(j in 1:ncol(phenotable)) {
			grid::grid.text(as.character(phenotable[i, j]),
					vp = grid::viewport(layout.pos.row = i + 1, layout.pos.col = j + 1))
		}
	}

	# Returning to the top viewport
	grid::upViewport()

	return(invisible(NULL))
}


#' Plotting haplotypes
#'
#' Details
#'
#' @param plotting_data To complete
#' @param difflist To complete
#' @param fontsize To complete
#'
#' @return To complete
#'
#' @export
#' @examples
#' NULL
grid.haplotypes <- function(plotting_data, difflist, fontsize = 8) {

	# Determining the number of columns as the maximum position to plot
	ncolumns <- max(do.call("rbind", plotting_data)$pos)

	# Dividing the viewport into rows
	# Haplotype rows are interleaved with rows used to show the differences
	#  between consecutive haplotypes, hence the number of rows being 2 * (number of haplotypes) - 1
	grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = length(plotting_data) * 2 - 1,
								     ncol = ncolumns)))

	# Plotting each haplotype in its row
	for(i in 1:length(plotting_data)) {
		hapdata <- plotting_data[[i]]

		for(j in 1:nrow(hapdata)) {
			grid::grid.text(hapdata[j, "nuc"],
					gp = grid::gpar(fontsize = fontsize,
							fontfamily = "mono",
							fontface = ifelse(hapdata[j, "log10p"] != 0, "bold", "plain"),
							col = hapdata[j, "color"]),
					vp = grid::viewport(layout.pos.row = i * 2 - 1,
							    layout.pos.col = hapdata[j, "pos"]))
		}
	}

	# Also plotting the positions where the nucleotides differ
	for(i in 1:length(difflist)) {
		if(!length(difflist[[i]])) next

		for(j in difflist[[i]]) {
			grid::grid.lines(x = 0.5, y = c(0, 1), vp = grid::viewport(layout.pos.row = i * 2, layout.pos.col = j))
		}
	}

	# Moving back into the top viewport
	grid::upViewport()

}

#' Find positions that differ between two sequences
#'
#' Details
#'
#' @param hapdata To complete
#'
#' @return To complete
#'
#' @export
#' @examples
#' NULL
nucdiff <- function(hapdata) {
	if(length(hapdata) <= 1) stop("nucdiff needs at least two haplotypes to compare")

	# Creating a list that will contain the positions that differ
	output <- list()

	for(i in 1:(length(hapdata) - 1)) {
		if(!all(hapdata[[i]]$pos == hapdata[[i + 1]]$pos)) stop("All positions must be shared between all haplotypes")

		output[[i]] <- which(hapdata[[i]]$nuc != hapdata[[i + 1]]$nuc & hapdata[[i]]$nuc != "-" & hapdata[[i + 1]]$nuc != "-")
	}

	output
}


#' Map a set of numeric values onto a color palette
#'
#' Details
#'
#' @param values To complete
#' @param max To complete
#' @param min To complete
#' @param pal To complete
#' @param n.colors To complete
#'
#' @return To complete
#'
#' @export
#' @examples
#' NULL
map_color <- function(values, max, min = 0, pal = "YlOrRd", n.colors = 9) {
	# Breaking the numeric values into a factor
	index <- cut(values, breaks = seq(min, max, length.out = n.colors))
	palette <- RColorBrewer::brewer.pal(n = n.colors, name = pal)
	output <- palette[as.integer(index)]
	# Setting the values at the minimum to black
	output[values <= min] <- "black"
	output
}

#' Adjust plotting position according to the presence of gaps
#'
#' Details
#'
#' @param hapdata To complete
#' @param alignment To complete
#'
#' @return To complete
#'
#' @export
#' @examples
#' NULL
adjust_gaps <- function(hapdata, alignment) {

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
#' Details
#'
#' @param x hapdata
#'
#' @return To complete
#'
#' @examples
#' NULL
fill_gaps <- function(hapdata) {
	# We need to get the set of positions covered from 1 to the maximum
	maxpos <- max(do.call("rbind", hapdata)$pos)

	# Then we loop over all the data.frames and fill the missing positions with dashes
	for(i in 1:length(hapdata)) {
		indices <- which(! 1:maxpos %in% hapdata[[i]]$pos)

		if(!length(indices)) next

		new_rows <- data.frame(pos = indices,
				       nuc = "-",
				       log10p = 0,
				       stringsAsFactors = FALSE)

		hapdata[[i]] <- rbind(hapdata[[i]], new_rows)
		hapdata[[i]] <- hapdata[[i]][order(hapdata[[i]]$pos), ]
	}

	hapdata
}

#' A function that reads the (sorted) k-mer p-values from the output of a k-mer GWAS analysis
#'
#' Details
#'
#' @param kmer_file To complete
#' @param max_kmers To complete
#' @param kmer_length To complete
#'
#' @return To complete
#'
#' @export
#' @examples
#' NULL
read_kmer_pvalues <- function(kmer_file, max_kmers, kmer_length) {

	# Reading the k-mer sequences and p-values from file, using only relevant columns
	output <- read.table(kmer_file,
			     colClasses = c("NULL", "character", rep("NULL", 6), "numeric"),
			     nrows = max_kmers)

	# Naming the two columns
	names(output) <- c("kmer_id", "pvalue")

	# As we assume that the p-value column is sorted, we need to check for this
	if(is.unsorted(output$pvalue)) stop("ERROR: read_kmer_pvalues expects sorted p-values in ", kmer_file)

	# Extracting the k-mer sequence from the k-mer ID
	output$kmer <- sub("[^ATGC]+", "", output$kmer_id)

	# Getting the reverse complement of the k-mer sequence
	output$kmer_reverse <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(output$kmer)))

	# Checking that all k-mer sequences have the expected length
	if(!all(nchar(output$kmer) == kmer_length) || !all(nchar(output$kmer_reverse) == 31)) {
		stop("Not all k-mers match expected length of ", kmer_length)
	}

	# Only keeping the necessary columns for returning the data.frame
	return(output[, c("kmer", "kmer_reverse", "pvalue")])
}

#' Reading the consensus sequences obtained at a locus for a set of samples
#'
#' Details
#'
#' @param input_fasta To complete
#'
#' @return To complete
#'
#' @export
#' @examples
#' NULL
read_consensus <- function(input_fasta) {
	fasta_lines <- readLines(input_fasta)

	# Extracting the names of the samples
	sample_names <- grep("^>", fasta_lines, value = TRUE)
	sample_names <- sub("^>", "", sample_names)

	# Extracting the sequences themselves from the fasta file
	sequences <- grep("^>", fasta_lines, value = TRUE, invert = TRUE)

	# Checking that the number of samples matches the number of sequences
	stopifnot(length(sample_names) == length(sequences))

	# The names of the samples are the names of the sequence vector
	names(sequences) <- sample_names

	return(sequences)
}

#' Extract the haplotypes found in a set of sequences and filter based on frequency
#'
#' Details
#'
#' @param sequences To complete
#' @param min_frequency To complete
#'
#' @return To complete
#'
#' @export
#' @examples
#' NULL
get_haplotypes <- function(sequences, min_frequency) {

	# First compute the number of times that a given sequence occurs in the dataset
	haplotypes <- table(sequences)

	# Keeping only haplotypes that occur at least min_frequency times
	haplotypes <- haplotypes[haplotypes >= min_frequency]

	# We return the sequences themselves instead of the table
	haplotypes <- names(haplotypes)
	haplotypes
}

#' A function that links haplotypes and the corresponding observed phenotypes
#'
#' Details
#'
#' @param sequences To complete
#' @param haplotypes To complete
#' @param phenotypes To complete
#' @param id_column To complete
#' @param phenotype_column To complete
#'
#' @return To complete
#'
#' @export
#' @examples
#' NULL
link_phenotypes <- function(sequences, haplotypes, phenotypes, id_column, phenotype_column) {

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

#' A function that matches a set of significant k-mers to positions on haplotypes
#'
#' Details
#'
#' @param haplotypes To complete
#' @param kmers To complete
#' @param kmer_length To complete
#'
#' @return To complete
#'
#' @export
#' @examples
#' NULL
match_kmers <- function(haplotypes, kmers, kmer_length) {

	# Creating a list that will have as many elements as there are haplotypes
	# This list will be output by the function
	overlap_list <- list()

	for(i in haplotypes) {
		overlap_list[[i]] <- kmers
		overlap_list[[i]]$fmatch   <- sapply(kmers$kmer, function(kmer) regexpr(kmer, i, fixed = TRUE))
		overlap_list[[i]]$rmatch   <- sapply(kmers$kmer_reverse, function(kmer) regexpr(kmer, i, fixed = TRUE)) 
		overlap_list[[i]]$matchpos <- pmax(overlap_list[[i]]$fmatch, overlap_list[[i]]$rmatch)
	}

	# Keeping only the k-mers for which there is at least one overlapping k-mer with a p-value
	# Also keeping only relevant columns
	overlap_list <- lapply(overlap_list, function(x) x[x$matchpos != -1, c("pvalue", "matchpos")])

	# Coercing to an IRanges object
	overlap_list <- lapply(overlap_list, function(x) {
				       IRanges::IRanges(start = x$matchpos, width = kmer_length, log10p = -log10(x$pvalue))
			     })

	overlap_list
}

#' Format haplotypes for plotting using k-mer overlap information
#'
#' Details
#'
#' @param haplotypes To complete
#' @param overlaps To complete
#'
#' @return To complete
#'
#' @export
#' @examples
#' NULL
format_haplotypes <- function(haplotypes, overlaps) {

	# Iterate over all haplotypes using lapply
	output <- lapply(haplotypes, function(x, kmer_overlaps) {

				 # Initializing the output data.frame with each nucleotide and its position
				 output_df <- data.frame(pos = 1:nchar(x),
							 nuc = strsplit(x, "")[[1]],
							 log10p = 0,
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

#' Align a set of haplotypes using mafft
#'
#' Details
#'
#' @param fasta_path To complete
#' @param haplotypes To complete
#' @param mafft_path To complete
#' @param mafft_options To complete
#'
#' @return To complete
#'
#' @export
#' @examples
#' NULL
mafft_align <- function(fasta_path, haplotypes, mafft_path, mafft_options) {

	# Opening the fasta file that will be used as input for mafft
	fasta <- file(fasta_path, open = "w+")
	on.exit(close(fasta))

	# Writing the haplotypes in fasta format to fasta_path
	for(i in 1:length(haplotypes)) {
		cat(paste0(">hap", i, "\n"), file = fasta)
		cat(haplotypes[i], "\n", file = fasta)
	}

	# Prepare the mafft command and launch it using the system interface
	# The output of mafft is read directly into an R variable
	mafft_command <- paste0(mafft_path, " ", mafft_options, " ", fasta_path)
	alignment <- system(mafft_command, intern = TRUE)

	# Formatting the alignment for output
	alignment <- strsplit(paste0(alignment, collapse = ""), split = ">hap[0-9]+")[[1]]
	alignment <- alignment[nchar(alignment) > 0]

	alignment
}

