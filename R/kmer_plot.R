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
	# Computing a table of the phenotypes observed per character
	phenotable <- table(haplotype_data$phenotype_character, haplotype_data$id)

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
#' @param gaps To complete
#'
#' @return To complete
#'
#' @export
#' @examples
#' NULL
adjust_gaps <- function(hapdata, gaps) {
	if(!length(gaps)) return(hapdata)

	for(i in 1:length(gaps)) {
		hapdata[hapdata$pos >= IRanges::start(gaps[i]), "pos"] <-
			hapdata[hapdata$pos >= IRanges::start(gaps[i]), "pos"] + IRanges::width(gaps[i])
	}

	return(hapdata)
}

#' Fill the deleted positions with dashes
#'
#' Details
#'
#' @param x To complete
#'
#' @return To complete
#'
#' @export
#' @examples
#' NULL
fill_gaps <- function(x) {
	# We need to get the set of positions covered from 1 to the maximum
	maxpos <- max(do.call("rbind", x)$pos)

	# Then we loop over all the data.frames and fill the missing positions with dashes
	for(i in 1:length(x)) {
		indices <- which(! 1:maxpos %in% x[[i]]$pos)
		if(!length(indices)) next
		new_rows <- data.frame(pos = indices,
				       nuc = "-",
				       log10p = 0,
				       color = "black",
				       stringsAsFactors = FALSE)
		x[[i]] <- rbind(x[[i]], new_rows)
		x[[i]] <- x[[i]][order(x[[i]]$pos), ]
	}

	x
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

