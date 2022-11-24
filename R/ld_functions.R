#' Check the validity of an LD matrix
#'
#' This is a convenience function for verifying the validity of an LD matrix
#' used by any function that takes such a matrix as input. It is not exported
#' by the package and therefore only meant to be used internally.
#'
#' @param ld_matrix A square numeric matrix with each cell indicating the pairwise
#'   linkage disequilibrium between two k-mers.
#'
#' @return A logical value indicating whether the input LD matrix is valid.
#' @examples
#' NULL
is_valid_ld <- function(ld_matrix) {

	if(!is.matrix(ld_matrix) || nrow(ld_matrix) != ncol(ld_matrix)) {
		stop("The input LD matrix should be a square matrix.")
	}

	if(!identical(rownames(ld_matrix), colnames(ld_matrix))) {
		stop("The row and column names of the LD matrix should be identical.")
	}

	if(!all(diag(ld_matrix) == 1)) {
		stop("All values along the diagonal of the LD matrix should be 1.")
	}

	if(any(is.na(ld_matrix))) {
		stop("NA values not supported within LD matrix.")
	}

	if(!all(ld_matrix >= 0 & ld_matrix <= 1)) {
		stop("All values of the LD matriix should lie between 0 and 1.")
	}

	return(TRUE)
}

#' Re-arrange the samples in an LD matrix according to user-specified criteria
#'
#' This function takes an LD matrix as input, a dataset of k-mer positions (with
#' associated p-values), and a sorting parameter, and returns the same LD matrix
#' with the order of the samples re-arranged according to the sorting parameter.
#'
#' @param ld_matrix A square numeric matrix with each cell indicating the pairwise
#'   linkage disequilibrium between two k-mers. The rows and columns of the matrix
#'   must be named after the k-mers.
#' @param kmer_data A GRanges object indicating the genomic location of significantly
#'   associated k-mers. This object must have its elements named after the k-mer
#'   sequence (for subsetting purposes) and have a metadata column named "log10p"
#'   for sorting based on p-values.
#' @param sort_param A parameter to use for sorting the LD matrix. At the moment
#'   two values are supported:
#'   \itemize{
#'     \item "position": k-mers are sorted according to their genomic position 
#'     \item "pvalue": k-mers are sorted from the most significant to least significant p-value
#'   }
#' @param positions A GRanges object used to restrict the LD matrix to the k-mers overlapping
#'   a given set of genomic positions. If \code{NULL} (the default), then no
#'   subsetting is performed.
#' 
#' @return An LD matrix that has been reordered according to the choice of sort_param,
#'   and possibly subset based on the positions parameter.
#'
#' @export
#' @examples
#' NULL
ld_sort <- function(ld_matrix, kmer_data, sort_param = c("position", "pvalue"), positions = NULL) {
	# Checking the validity of the input LD matrix.
	stopifnot(is_valid_ld(ld_matrix))

	# Checking the validity of the k-mer data input
	if(! "log10p" %in% names(S4Vectors::mcols(kmer_data))) {
		stop("Error in ld_sort: Metadata column 'log10p' not found in kmer_data")
	}

	# First we need to check that all the k-mers in the LD matrix are in kmer_data
	if(!all(rownames(ld_matrix) %in% names(kmer_data))) {
		stop("Error in ld_sort: Not all k-mers in ld_matrix are represented in kmer_data.\n",
		     "Ensure that kmer_data elements are named after k-mer sequence.")
	}

	# Subsetting the input kmer_data using the positions parameter (if not NULL)
	if(!is.null(positions)) {
		stopifnot(inherits(positions, "Granges"))
		kmer_data <- IRanges::subsetByOverlaps(kmer_data, positions, ignore.strand = TRUE)
	}

	# Next we re-order the kmer_data object based on the parameter passed
	if(!sort_param[1] %in% c("position", "pvalue")) {
		stop("Error in ld_sort: sort_param must be one of 'position' or 'pvalue'")
	}

	if(sort_param[1] == "position") {
		kmer_data <- sort(kmer_data, ignore.strand = TRUE)
	} else if(sort_param[1] == "pvalue") {
		kmer_data <- kmer_data[order(kmer_data$log10p, decreasing = TRUE)]
	}

	# We only keep the k-mers that are in the input matrix
	kmer_data <- kmer_data[names(kmer_data) %in% rownames(ld_matrix)]

	# Now we proceed to the reordering
	ld_matrix <- ld_matrix[, names(kmer_data)]
	ld_matrix <- ld_matrix[names(kmer_data), ]

	return(ld_matrix)
}

#' Subsample from a set of significant k-mers prior to computing LD
#'
#' Computing the pairwise LD between all k-mers in a set is computationally
#' expensive as the number of comparisons made is the square of the number
#' of input k-mers. Therefore, this function allows for subsampling from
#' a set of input k-mers such that the most significant k-mers are represented
#' while sampling from the whole genome and not just from the most significant
#' regions. This is acheived from selecting a fixed number of top significant
#' k-mers (the npvalue parameter) and sampling from the remaining k-mers
#' with a probability that is inversely proportional to the number of k-mers
#' found on a given sequence. This ensures that k-mers are sampled from
#' various regions of the genome.
#'
#' @param kmer_data A GRanges object containing the genomic positions of 
#'   significantly associated k-mers. Must containg a "log10p" column
#'   to sample k-mers based on their p-values, and a "kmer_canon" column
#'   containing the canonized sequence of the k-mer matching at that position.
#' @param nkmers A numeric indicating the total number of k-mers to subsample.
#' @param npvalue A numeric indicating the number of k-mers to sample from
#'   the top p-values. The npvalue most associated k-mers will be sampled prior
#'   to random sampling of the remaining k-mers.
#'
#' @return A subset of the input kmer_data object containing the only the
#'   subsampled k-mers.
#'
#' @export
#' @examples
#' NULL
subsample_kmers <- function(kmer_data, nkmers, npvalue) {

	# Checking the validity of the inputs
	stopifnot(inherits(kmer_data, "GRanges"))

	if(!all(c("log10p", "kmer_canon") %in% names(S4Vectors::mcols(kmer_data)))) {
		stop("Error in subsample_kmers: the input kmer_data object must include 'log10p' and 'kmer_canon' metadata columns.")
	}

	if(nkmers > length(kmer_data)) {
		stop("Error in subsample_kmers: the number of k-mers to sample (", nkmers,
		     ") is greater than the number of k-mers in kmer_data (", length(kmer_data), ")")
	}

	if(npvalue >= nkmers) {
		stop("Error in subsample_kmers: npvalue must be strictly smaller than nkmers.")
	}

	# Computing the number of k-mers found on each sequence
	seqnames_table <- table(GenomicRanges::seqnames(kmer_data))
	kmer_data$seq_count <- seqnames_table[as.character(GenomicRanges::seqnames(kmer_data))]

	# Subsetting a number of k-mers based on their p-values
	kmer_data <- kmer_data[order(kmer_data$log10p, decreasing = TRUE)]
	pvalue_kmers <- kmer_data$kmer_canon[1:npvalue]

	# We need to pick the remaining k-mers from the ones that have not been chosen yet
	remaining_kmers <- kmer_data[!kmer_data$kmer_canon %in% pvalue_kmers]
	nothers <- nkmers - npvalue

	# We sample from the remaining k-mers with probability that is inversely proportional
	# to the number of k-mers matching a given sequence. This way we get to see k-mers that
	# are located on chromosomes other than the one on which the main signal is found
	other_kmers <- sample(remaining_kmers$kmer_canon, nothers, replace = FALSE, prob = 1 / remaining_kmers$seq_count)

	# The output k-mers are a combination of the ones from the p-values and the randomly selected ones
	kmer_data[kmer_data$kmer_canon %in% c(pvalue_kmers, other_kmers)]
}

#' Select colors for a discrete scale as in ggplot2
#'
#' This code is taken from https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#' It is meant to be used internally by the package and is therefore not exported.
#'
#' @param n A numeric value. The number of colors to generate.
#'
#' @return A vector of colors to use for plotting.
#'
#' @examples
#' NULL
gg_color_hue <- function(n) {
	hues = seq(15, 375, length = n + 1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Plot an LD matrix using grid functions
#'
#' This function uses grid functionality to plot an LD matrix.
#' LD values are represented using different shades from yellow
#' to red as implemented by the RColorBrewer "YlOrRd" scale.
#' The chromosome or scaffold to which a k-mer belong is represented
#' by colored rectangles plotted in the left outer margin of the plot.
#' Similarly, the -log10(p-values) are represented using the "Blues"
#' color scale of the package RColorBrewer in the bottom outer
#' margin of the plot. The y-axis of the generated plot is inverted
#' with respect to the input matrix, i.e., the first row of the matrix
#' is found at the very bottom of the plot and the last row of the matrix
#' is found at the very top of the plot.
#'
#' @param ld_matrix A matrix of pairwise linkage disequilibrium values observed
#'   among a set of k-mers.
#' @param kmer_positions A Granges object storing the genomic positions and
#'   p-values associated with a set of significant k-mers. Must contain
#'   a 'log10p' columns containing the -log10(p-value) associated with each k-mer.
#'   All k-mers represented in the ld_matrix object must also be presented within
#'   kmer_positions.
#' @param top_legend A logical. Whether or not to plot the color scale for the
#'   chromosomes at the top of the plot. Defaults to TRUE.
#' @param ylabels A logical. If TRUE, labels are printed on the y-axis at the mid-position
#'   for each reference sequence to identify them in the plot. Setting this option to
#'   TRUE only makes sense when k-mers are sorted by genomic position, therefore it is
#'   FALSE by default.
#' @param ylabel_pattern A regular expression to use to keep only the labels matching
#'   that pattern. This option is only relevant if ylabels is set to TRUE.
#' @param fontsize A numeric. The font size to use for the y-labels and color scales.
#'
#' @return \code{NULL}, invisibly. This function is invoked for its plotting side effect.
#'
#' @export
#' @examples
#' NULL
ld_plot <- function(ld_matrix, kmer_positions, top_legend = TRUE, ylabels = FALSE,
		    ylabel_pattern = NULL, fontsize = 8) {

	# Checking the validity of the ld_matrix input
	stopifnot(is_valid_ld(ld_matrix))

	# Checking the validity of the kmer_positions object
	if(! inherits(kmer_positions, "GRanges") || ! "log10p" %in% names(S4Vectors::mcols(kmer_positions))) {
		stop("Error in ld_plot: kmer_positions must be a GRanges object with a 'log10p' column.")
	}

	if(!all(rownames(ld_matrix) %in% names(kmer_positions))) {
		stop("Error in ld_plot: Not all k-mers in ld_matrix are represented in kmer_positions\n",
		     "Ensure that kmer_positions elements are named after k-mer sequence.")
	}

	# Map the numeric values onto a color scale
	color_scale <- map_color(values = ld_matrix, pal = "YlOrRd",
				 n_colors = 9, scale_extend = 0, n_colors_remove = 0)
	color_values <- matrix(color_scale$mapped_colors, nrow = nrow(ld_matrix))

	# Also preparing a vector of colors to show the chromosome that a k-mer belongs to
	kmer_chrom <- as.character(GenomicRanges::seqnames(kmer_positions[rownames(ld_matrix)]))
	n_chrom <- length(unique(kmer_chrom))

	legend_text  <- unique(kmer_chrom)
	legend_color <- gg_color_hue(n_chrom)[as.factor(legend_text)]
	names(legend_color) <- legend_text
	legend_color <- legend_color[order(names(legend_color))]

	kmer_color <- legend_color[kmer_chrom]

	# Doing the same thing by mapping the p-values on a blue color-scale on the x-axis
	kmer_pvalues  <- kmer_positions[colnames(ld_matrix)]$log10p
	pvalue_scale  <- map_color(values = kmer_pvalues, pal = "Blues",
				   n_colors = 9, scale_extend = 0, n_colors_remove = 0)
	pvalue_colors <- pvalue_scale$mapped_colors

	# Creating a viewport layout suitable for out plotting purposes
	grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 3,
								     ncol = 3,
								     widths = grid::unit(c(1, 0.75, 1), c("null", "npc", "null")),
								     heights = grid::unit(c(1, 0.75, 1), c("null", "npc", "null")))))

	# Push a viewport with suitable coordinates for positioning the rectangles to plot
	grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 2,
					  xscale = c(1, ncol(ld_matrix) + 1),
					  yscale = c(1, nrow(ld_matrix) + 1)))

	# Formatting the data for plotting
	plotting_data <- expand.grid(i = 1:nrow(ld_matrix), j = 1:ncol(ld_matrix))
	plotting_data$color <- color_values[as.matrix(plotting_data)]

	# Plotting using rectangles in the appropriate viewports
	grid::grid.rect(x = plotting_data$j,
			y = plotting_data$i,
			width = 1,
			height = 1,
			just = c(0, 0),
			default.units = "native",
			gp = grid::gpar(col = "transparent", fill = plotting_data$color))

	# Adding colored lines to the left side of the plot indicating the chromosome on which a k-mer is found
	# We go to the left viewport for that
	grid::upViewport()
	grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1,
					  yscale = c(1, nrow(ld_matrix) + 1)))

	grid::grid.rect(y = 1:length(kmer_color),
			x = grid::unit(0.95, "npc"),
			height = 1,
			width = grid::unit(0.5, "npc"),
			default.units = "native",
			just = c(1,0),
			gp = grid::gpar(col = "transparent", fill = kmer_color))

	# Optionally adding y-axis labels for the reference sequence
	if(ylabels) {
		for(i in legend_text) {
			if(!is.null(ylabel_pattern) && !grepl(ylabel_pattern, i)) next
			ypos <- which(kmer_chrom == i)
			median_pos <- median(ypos)
			min_pos <- min(ypos)

			# Adding the label
			grid::grid.text(i,
					x = grid::unit(0.38, "npc"),
					y = grid::unit(median_pos, "native"),
					just = "right",
					gp = grid::gpar(fontsize = fontsize))

			# Adding a horizontal line to separate neighouring chromosomes
			grid::grid.lines(x = grid::unit(c(0.45, 0.95), "npc"),
					 y = grid::unit(min_pos, "native"),
					 gp = grid::gpar(lwd = 0.5))
		}
	}


	# Adding colored lines to below the plot indicating the -log10(p-value) of the k-mer
	# We go to the botton viewport for that
	grid::upViewport()
	grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 2,
					  xscale = c(1, ncol(ld_matrix) + 1)))
	grid::grid.rect(y = grid::unit(0.95, "npc"),
			x = 1:length(pvalue_colors),
			height = grid::unit(0.5, "npc"),
			width = 1,
			default.units = "native",
			just = c(0,1),
			gp = grid::gpar(col = "transparent", fill = pvalue_colors))

	# Optionally going to the top viewport to plot the color legend for the chromosomes
	if(top_legend) {
		grid::upViewport()
		grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2,
						  layout = grid::grid.layout(nrow = 1, ncol = length(legend_color))))

		for(i in 1:length(legend_color)) {
			grid::grid.rect(x = 0, width = 0.3, height = 0.3, just = "left",
					gp = grid::gpar(col = "transparent", fill = legend_color[i]),
					vp = grid::viewport(layout.pos.col = i))
			grid::grid.text(names(legend_color)[i], x = 0.4, just = "left",
					gp = grid::gpar(fontsize = fontsize),
					vp = grid::viewport(layout.pos.col = i))
		}
	}

	grid::upViewport()

	# Plotting the color scale for the LD values in the right viewport (top half)
	grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 3))
	grid::pushViewport(grid::viewport(x = 0, y = 0.5, height = 0.5, width = 0.5, just = c(0, 0)))

	grid.colorscale(breaks = color_scale$breaks,
			base_palette = color_scale$base_palette,
			label_text = expression(LD (r^2)),
			round_digits = 2,
			direction = "vertical",
			fontsize = fontsize)

	grid::upViewport()

	# Doing the same thing for the p-value scale in the lower half of the right viewport
	grid::pushViewport(grid::viewport(x = 0, y = 0, height = 0.5, width = 0.5, just = c(0, 0)))

	grid.colorscale(breaks = pvalue_scale$breaks,
			base_palette = pvalue_scale$base_palette,
			label_text = expression(-log[10](italic(p))),
			round_digits = 1,
			direction = "vertical",
			fontsize = fontsize)

	# Going back up to the layout viewport
	grid::upViewport()

	# And then to the top viewport
	grid::upViewport()

	invisible(NULL)
}

#' Compute the pairwise LD between k-mers based on their presence/absence
#'
#' This function has been obsoleted by a much faster C implementation
#'
#' @param pav A matrix with presence/absence of a set of k-mers
#' @param kmers A character vector of the k-mer sequences
#'
#' @return A numeric matrix with the pairwise LD values computed between
#'   all k-mers.
#' @examples
#' NULL
kmer_ld <- function(pav, kmers) {
	# Creating an output matrix of N x N
	output_matrix <- matrix(nrow = nrow(pav), ncol = nrow(pav))

	# The diagonal is 1 by definition
	diag(output_matrix) <- 1

	# Then we loop over the kmers
	for(i in 1:(nrow(pav) - 1)) {
		for(j in (i + 1):nrow(pav)) {
			if(i == j) next

			# Computing the LD for this pair of k-mers
			pxy1 <-  sum(pav[i, ] == 1 & pav[j, ] == 1) / ncol(pav)
			px1 <- sum(pav[i, ] == 1) / ncol(pav)
			py1 <- sum(pav[j, ] == 1) / ncol(pav)
			px0 <- 1 - px1
			py0 <- 1 - py1

			d <- pxy1 - px1 * py1
			denominator <- px1 * py1 * px0 * py0

			r2 <- d^2 / denominator

			output_matrix[i, j] <- r2
		}
	}

	output_matrix[lower.tri(output_matrix)] <- t(output_matrix)[lower.tri(output_matrix)]
	rownames(output_matrix) <- colnames(output_matrix) <- kmers
	output_matrix
} 


#' Greedy clustering of an LD matrix
#'
#' This function has been obsoleted by a much faster C implementation
#'
#' @param ld_matrix A matrix of pairwise LD values between a set of k-mers
#'
#' @return A numeric matrix with the pairwise LD values computed between all
#'   k-mers, clustered such that k-mers with high LD are grouped together.
#'
#' @examples
#' NULL
cluster_ld <- function(ld_matrix) {
	to_cluster <- rownames(ld_matrix)[-1]
	clustered <- rownames(ld_matrix)[1]

	while(length(to_cluster) > 0) {
		current_sample <- clustered[length(clustered)]
		closest_sample <- which.max(ld_matrix[current_sample, to_cluster])
		clustered <- c(clustered, to_cluster[closest_sample])
		to_cluster <- to_cluster[-closest_sample]
	}

	# Creating a matrix with the same dimensions as the input one
	output_matrix <- matrix(nrow = nrow(ld_matrix), ncol = ncol(ld_matrix))	
	rownames(output_matrix) <- colnames(output_matrix) <- clustered

	# Filling the values in the new matrix using the ones from the input matrix
	for(i in rownames(output_matrix)) {
		for(j in colnames(output_matrix)) {
			output_matrix[i, j] <- ld_matrix[i, j]
		}
	}

	output_matrix
}

