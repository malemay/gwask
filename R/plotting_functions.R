#' manhattan_plot
#' 
#' @param formatted_data A GRanges object formatted with either \code{\link{format_kmer_gwas}}
#'   or \code{\link{format_gapit_gwas}} depending on the GWAS type. These functions
#'   return GRanges containing standard columns and seqinfo in addition to GWAS type-specific
#'   columns, therefore allowing a common interface for all GWAS types.
#' @param gwas_type A character of length one. Only two options are supported at
#'   this time, either \code{"kmer"} for k-mer GWAS as implemented by Voichek
#'   and Weigel (2020) or \code{"gapit"} for the output of GAPIT.
#' @param signals A data.frame with the coordinates of GWAS signals to 
#'   show on the plot. If \code{NULL} (the default) or \code{nrow(signals)} < 0,
#'   then no signals are plotted. If supplied, the data.frame must minimally
#'   contain the following columns:
#'   \itemize{
#'     \item  first_chrom: the chromosome where the leftmost position is found
#'     \item  first_pos: the coordinate on first_chrom where the leftmost significant marker occurs
#'     \item  last_chrom: the chromosome where the rightmost position is found
#'     \item  last_pos: the coordinate on last_chrom where the rightmost significant marker occurs
#'     \item  max_chrom: the chromosome where the marker with the most significant p-value occurds
#'     \item  max_pos: the coordinate on max_chrom where the most significant marker occurs
#'   }
#'   Supplying this data.frame in itself only results in a vertical line being
#'   drawn at the location of each signal. If, in addition, \code{extent = TRUE},
#'   then the extent of the signal (from the first to the last position along
#'   the chromosome) is also plotted as a rectangle.
#' @param extent Logical. Whether or not to plot the extent of the signals
#'   as a rectangle. Can only be \code{TRUE} if signals is not \code{NULL}.
#' @param single_panel Logical. Whether all chromosomes should be plotted in
#'   a single panel (the default) or \code{\link[ggplot]{facet_wrap}} should
#'   be used to plot each chromosome in a separate panel.
#' @param threshold Numeric. A single value indicating where the -log10(p-value)
#'   threshold line should be plotted.
#' @param min_log10p Numeric. The minimum -log10(p-value) to display on the plot
#'   can be used to filter out some non-significant markers for plotting speedup.
#' @param plot_mapq Logical. Whether to plot the mapping quality (MAPQ) as the
#'   color aesthetic (when provided as k-mers data).
#'
#' @return A ggplot object describing a Manhattan plot with the specified parameters.
#'
#' @export
#' @examples
#' NULL
manhattan_plot <- function(formatted_data, gwas_type,
			   signals = NULL, extent = FALSE, single_panel = TRUE,
			   threshold = NULL, min_log10p = 0, plot_mapq = FALSE) {

	# Sanity checks -------------------------------------

	# Checking that formatted_data has the right format
	if(!"log10p" %in% names(mcols(formatted_data))) {
		stop("The formatted_data object does not satisfy format requirements.")
	}

	# Checking that gwas_type is supported
	if(!gwas_type %in% c("kmer", "gapit")) {
		stop("gwas_type ", gwas_type, " is not supported.")
	}

	# Checking that the MAPQ column is suppled if type == kmer
	if(gwas_type == "kmer" && ! "mapq" %in% names(mcols(formatted_data))) {
		stop("mapq column must be present in formatted_data for gwas_type = \"kmer\"")
	}

	# Checking that signals has the right format (if supplied)
	if(!is.null(signals) && !all(c("first_chrom", "first_pos", "last_chrom", "last_pos", "max_chrom", "max_pos") %in% names(signals))) {
		stop("The signals data.frame does not satisfy format requirements.")
	}

	# Checking that extent is TRUE only if signals is supplied
	if(is.null(signals) && extent) {
		stop("extent can only be TRUE if signals is not NULL")	
	}

	# Check that mapq is only TRUE if the MAPQ column actually exists
	if(plot_mapq && ! "MAPQ" %in% names(mcols(formatted_data))) {
		warning("plot_mapq = TRUE but \"MAPQ\" is not among the columns. Setting plot_mapq to FALSE")
		plot_mapq <- FALSE
	}

	# End of sanity checks ------------------------------

	# Using the seqinfo of the formatted_data input to compute various information
	chrom_lengths <- GenomeInfoDb::seqlengths(formatted_data)

	chrom_start <- cumsum(chrom_lengths)
	chrom_start <- c(0, chrom_start[-(length(chrom_start))])

	label_pos <- chrom_start + 0.5 * chrom_lengths

	# All vectors have the names of the sequences to allow lookup of the info by name
	names(chrom_start) <- names(chrom_lengths)  <- names(label_pos) <- GenomeInfoDb::seqlevels(formatted_data)

	# Computing the absolute position on the genome
	formatted_data$manhattan_rpos <- GenomicRanges::start(formatted_data) + chrom_start[as.character(GenomicRanges::seqnames(formatted_data))]

	# Computing whether this chromsome is even or odd (for plotting purposes)
	even_chromosomes <- GenomeInfoDb::seqlevels(formatted_data)[c(FALSE, TRUE)]
	formatted_data$manhattan_even <- as.character(GenomicRanges::seqnames(formatted_data)) %in% even_chromosomes

	# Issuing a warning until this functionality is restored
	warning("signals functionality in manhattan_plot() not functional")

	# Formatting the signals data for plotting (the signals format should be redesigned)
	if(!is.null(signals) && nrow(signals)) {
		signals$first_plot <- fai_info$start[signals$first_chrom] + signals$first_pos
		signals$last_plot <- fai_info$start[signals$last_chrom] + signals$last_pos
		signals$max_plot <- fai_info$start[signals$max_chrom] + signals$max_pos
		# For compatibility with facet_wrap when extent = TRUE
		signals$manhattan_chrom <- signals$first_chrom
	}

	# Removing the associations that are less significant than min_log10p
	if(min_log10p > 0) {
		formatted_data <- formatted_data[formatted_data$log10p >= min_log10p]
	}

	# If single_panel is TRUE, then all the data is plotted along a single axis
	if(single_panel) {
		output <- ggplot2::ggplot(as.data.frame(formatted_data)) +
			ggplot2::geom_point(mapping = ggplot2::aes_string(x = "manhattan_rpos",
									  y = "log10p",
									  color = ifelse(plot_mapq, "MAPQ", "manhattan_even"))) +
				ggplot2::scale_x_continuous(name = "Chromosome",
							    breaks = label_pos,
							    labels = names(label_pos),
							    limits = c(0, sum(chrom_lengths))) +
				ggplot2::ylab("-log10(p-value)") +
				ggplot2::theme_bw() +
				ggplot2::theme(panel.grid.minor = ggplot2::element_blank())


	} else {
	# If single_panel is not TRUE, there is one panel per chromosome
		output <- ggplot2::ggplot(formatted_data) +
			ggplot2::geom_point(mapping = ggplot2::aes_string(x = "start",
									  y = "log10p")) +
			ggplot2::facet_wrap(~seqnames, ncol = 5) +
			ggplot2::scale_x_continuous(name = "Position (bp)") +
			ggplot2::ylab("-log10(p-value)") +
			ggplot2::theme_bw() +
			ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
	}

	if(!is.null(signals) && nrow(signals)) {
		output <- output + ggplot2::geom_vline(data = signals, aes(xintercept = max_plot), linetype = 2)

		if(extent) {
			output <- output + ggplot2::geom_rect(data = signals,
							      ymin = 0, ymax = max(formatted_data$log10p),
							      ggplot2::aes(xmin = first_pos, xmax = last_pos),
							      color = "red", alpha = 0.4)
		}
	}

	# Unless plot_mapq is TRUE, we do not need the color scale
	if(!plot_mapq) output <- output + ggplot2::guides(color = "none")

	# Also adding the threshold if not NULL
	if(!is.null(threshold)) output <- output + ggplot2::geom_hline(yintercept = threshold, linetype = 2)

	return(output)
}

#' Generate a grob representing a transcript using grid functions
#'
#' This function generates a grob of a transcript model using grid functions.
#' 
#' The function assumes that that the plotting will occur in a viewport
#' with x-axis native coordinates already suitable for plotting that
#' transcript. y-axis coordinates do not matter and "npc" coordinates
#' are used for centering the transcript models vertically in the
#' viewport.
#'
#' The function will be typically called by the wrapper function
#' \code{\link{plot_transcripts}} that will take care of organizing
#' the viewports and coordinate scales for plotting all the transcripts
#' of a given gene model.
#'
#' @param tx_gene A GRanges object representing the coordinates of
#'   of a single gene to which the transcript belongs.
#' @param tx_exons A GRanges object representing the exons belonging
#'   to a single transcript.
#' @param tx_name A GRanges object representing the coding sequences
#'   belonging to a single transcript.
#'
#' @return A gList grob object representing a single transcript.
#'
#' @export
#' @examples
#' NULL
transcriptGrob <- function(tx_gene, tx_exons, tx_cds) {

	# Drawing the line that goes from start to end of the gene
	tx_line <- grid::linesGrob(x = grid::unit(c(GenomicRanges::start(tx_gene), GenomicRanges::end(tx_gene)), "native"),
				   y = grid::unit(0.5, "npc"),
				   name = "tx_line")

	# Drawing the rectangles that correspond to exons
	tx_exons <- grid::rectGrob(x = grid::unit(GenomicRanges::start(tx_exons), "native"),
				   y = grid::unit(0.5, "npc"),
				   width = grid::unit(GenomicRanges::width(tx_exons), "native"),
				   height = grid::unit(0.8, "npc"),
				   hjust = 0,
				   gp = grid::gpar(fill = "white"),
				   name = "tx_exons")

	# Doing the same thing for the coding sequences by using a shaded color
	tx_cds <- grid::rectGrob(x = grid::unit(GenomicRanges::start(tx_cds), "native"),
				 y = grid::unit(0.5, "npc"),
				 width = grid::unit(GenomicRanges::width(tx_cds), "native"),
				 height = grid::unit(0.8, "npc"),
				 hjust = 0,
				 gp = grid::gpar(fill = ifelse(as.character(GenomicRanges::strand(tx_cds)) == "+",
							       "skyblue", 
							       "orange")),
				 name = "tx_cds")

	return(gList(tx_line, tx_exons, tx_cds))
}

#' Plot all the possible transcripts for genes in a genomic region using grid functions
#'
#' This function calls plot_transcript on all individual transcripts
#' of the genes located in a given genomic region and plots each separately
#' in a row of a viewport grid.
#'
#' This function requires as input the full set of genes, transcripts,
#' exons and coding sequences for the reference genome being used. It 
#' will check that transcripts exist for all genes provided. This
#' behaviour may be modified in the future depending on use cases.
#'
#' It is also an error in the current implementation to try and plot
#' over an interval that does not overlap any genes.
#'
#' The function itself does not plot any objects but instead arranges
#' the viewport layout and moving between viewports, and calls the
#' \code{\link{plot_transcript}} to plot the transcripts one by one.
#' The latter is the function that should be modified for changing
#' the way that plotted transcripts look.
#'
#' The number of rows in the plotting viewport is the number of trancripts
#' for the gene with the most transcripts over the interval. This is done
#' so that enough space is provided to plot all transcripts for a given
#' gene below each other.
#'
#' @param genes A GRanges object containing the known genes for the
#'   reference genome being studied, as obtained from processing
#'   a TxDb object of the GenomicFeatures package.
#' @param transcripts A GRangesList object containing the possible
#'   transcripts for the set of genes in the reference genome being used.
#'   Each element of the GRangesList object is a GRanges object containing
#'   the known transcripts for a particular gene, with the name of the
#'   list element being the name of the gene.
#' @param exons A GRangesList object containing the documented exons
#'   in the reference genome being used. Each element of the list contains
#'   the exons for a particular transcript, with the name of the list
#'   element being the name of the transcript.
#' @param cds A GRangesList object containing the documented coding sequences
#'   in the reference genome being used. Each element of the list contains
#'   the coding sequences for a particular transcript, with the name of the list
#'   element being the name of the transcript.
#' @param xscale A GRanges object used to set the limits of the x-axis
#'   and to subset the genes object such that only the ones located in the 
#'   interval are plotted.
#' @param transcript_margins A numeric of length 4 or that can be recycled to that
#'   length. The margins used by the function \code{\link[grid]{plotViewport}}
#'   to set the margins around the plotting region. By default it is set
#'   to c(0, 4.1, 0, 2.1).
#'
#' @return NULL, invisibly. The function is called for its plotting
#'   side-effect.
#'
#' @export
#' @examples
#' NULL
plot_transcripts <- function(genes, transcripts, exons, cds, xscale,
			     transcript_margins = c(0, 4.1, 0, 2.1)) {

	# Checking that the inputs are of the right type
	stopifnot(inherits(genes, "GRanges"))
	stopifnot(inherits(transcripts, "GRangesList"))
	stopifnot(inherits(exons, "GRangesList"))
	stopifnot(inherits(cds, "GRangesList"))

	# Doing a few more sanity checks on the inputs
	stopifnot(length(xscale) == 1)
	stopifnot(all(names(genes) %in% names(transcripts)))
	stopifnot(sum(lengths(transcripts)) == length(exons))
	stopifnot(sum(lengths(transcripts)) == length(cds))

	# Extracting the genes to plot from the xscale GRanges object
	genes <- IRanges::subsetByOverlaps(genes, xscale)
	if(!length(genes)) stop("No genes in the selected interval")

	# Extracting the transcripts related to those genes
	transcripts <- transcripts[names(genes)]

	# The number of rows for the plotting viewports depends on the gene with the most transcripts
	nrows <- max(lengths(transcripts))

	# Pushing a viewport layout with as many rows as there are transcripts for the gene with the highest number of transcripts
	grid::pushViewport(grid::plotViewport(margins = transcript_margins, layout = grid::grid.layout(nrow = nrows)))

	# Iterating over the genes
	for(gene_index in 1:length(genes)) {

		# Getting the names of the transcripts of that gene
		tnames <- transcripts[[ names(genes)[gene_index] ]]$tx_name

		# Looping over all transcripts to plot them in separate viewports
		for(i in 1:length(tnames)) {
			tx_name <- tnames[i]
			grid::pushViewport(grid::plotViewport(layout.pos.row = i,
							      xscale = c(GenomicRanges::start(xscale),
									 GenomicRanges::end(xscale))))

			# Drawing a transcriptGrob for each transcript
			grid::grid.draw(transcriptGrob(genes[gene_index], exons[[tx_name]], cds[[tx_name]]))

			# Moving back to the parent viewport
			grid::upViewport()
		}
	}

	# Going back to the viewport that we started from
	grid::upViewport()

	return(invisible(NULL))
}

#' A grob representing p-values of a GWAS analysis at a locus to plot using grid functions
#'
#' This function takes a data.frame of formatted GWAS results and returns
#' a gTree object representing a p-value plot in a given interval.
#'
#' @param gwas_results A GRanges object or data.frame of GWAS results.
#'   If a data.frame, must include the columns manhattan_chrom and
#'   manhattan_cpos for coercion to a GRanges. The metadata column
#'   log10p must be present as it will be used to plot the
#'   p-values on the y-axis.
#' @param interval A GRanges object used to subset the plotting to a given
#'   region.
#' @param feature A GRanges object indicating the position of a feature of
#'  interest to mark with vertical dotted lines. Both the start and end of
#'  the feature will be indicated with a line. If NULL (default), then no
#'  lines are drawn.
#' @param pvalue_margins A numeric of length 4 or that can be recycled to that
#'   length. The margins used by the function \code{\link[grid]{plotViewport}}
#'   to set the margins around the plotting region. By default it is set
#'   to c(5.1, 4.1, 4.1, 2.1).
#' @param yexpand A numeric of length 2 representing expansion factors
#'   of y-scale limits. The first value is the expansion factor to
#'   the top while the second value is the expansion factor to the bottom.
#'   Values of 0 (the default) represent no expansion relative to the
#'   values in x-scale.
#'
#' @return A gTree object that can be plotted with grid.draw() to display the p-values.
#'
#' @export
#' @examples
#' NULL
pvalueGrob <- function(gwas_results, interval, feature = NULL,
			pvalue_margins  = c(5.1, 4.1 , 4.1, 2.1),
			yexpand = c(0, 0)) {

	if(!inherits(gwas_results, "GRanges")) {
		gwas_results <- GenomicRanges::makeGRangesFromDataFrame(gwas_results,
									keep.extra.columns = TRUE,
									ignore.strand = TRUE,
									seqnames.field = "manhattan_chrom",
									start.field = "manhattan_cpos",
									end.field = "manhattan_cpos")
	}

	# We keep only the part of the gwas_results that overlaps the interval
	gwas_results <- IRanges::subsetByOverlaps(gwas_results, interval)

	# Setting the yscale depending on whether GWAS data is provided or not
	if(length(gwas_results)) {
		# Setting the yscale interval
		stopifnot(length(yexpand) == 2)
		yrange <- max(gwas_results$log10p) - min(gwas_results$log10p)
		yscale <- c(min(gwas_results$log10p) - yexpand[1] * yrange,
			    max(gwas_results$log10p) + yexpand[2] * yrange)
	} else  {
		yscale <- c(0, 1)
	}

	# Creating a viewport with appropriate scales
	pvalue_viewport <- grid::plotViewport(margins = pvalue_margins,
					      xscale = c(GenomicRanges::start(interval),
							 GenomicRanges::end(interval)),
					      yscale = yscale)

	# Creating the output gTree object to which grobs will be added
	output_gtree <- grid::gTree(name = "pvalueGrob", vp = pvalue_viewport, cl = "pvalue")

	# We plot a box around the viewport and an x-axis
	output_gtree <- grid::addGrob(output_gtree, grid::rectGrob())
	output_gtree <- grid::addGrob(output_gtree, grid::xaxisGrob())
	output_gtree <- grid::addGrob(output_gtree, grid::textGrob("Position along reference (bp)", y = grid::unit(-3, "lines")))

	# Adding vertical lines with the location of the feature if provided
	if(!is.null(feature)) {
		output_gtree <- grid::addGrob(output_gtree,
					      grid::linesGrob(x = GenomicRanges::start(feature),
							      y = grid::unit(c(0, 1), "npc"),
							      default.units = "native",
							      gp = grid::gpar(lty = 2)))

		output_gtree <- grid::addGrob(output_gtree,
					      grid::linesGrob(x = GenomicRanges::end(feature),
							      y = grid::unit(c(0, 1), "npc"),
							      default.units = "native",
							      gp = grid::gpar(lty = 2)))
	}

	# An empty GWAS dataset over the range needs to be handled in a special way
	if(!length(gwas_results)) {
		warning("No GWAS data in interval")

		# Plotting some very basic features
		output_gtree <- grid::addGrob(output_gtree, grid::textGrob("No GWAS data in selected range"))

		# Returning from the function because the rest of the function should not be processed
		return(output_gtree)
	}

	# Then we can plot the data points
	output_gtree <- grid::addGrob(output_gtree,
				      grid::pointsGrob(x = grid::unit(GenomicRanges::start(gwas_results), "native"),
						       y = grid::unit(gwas_results$log10p, "native")))

	# Adding a y-axis
	output_gtree <- grid::addGrob(output_gtree, grid::yaxisGrob())

	# Adding the axis labels
	output_gtree <- grid::addGrob(output_gtree,
				      grid::textGrob("-log10(p-value)", x = grid::unit(-3, "lines"), rot = 90))

	return(output_gtree)
}

#' Plotting p-values along transcript model using grid functions
#'
#' This function plots the p-values observed along the genome in
#' the bottom panel and a model of the transcripts for genes in the
#' region in the top panel.
#'
#' The margins of the transcript and p-value plots are set separately
#' but the left and right margins of each must align in order for the
#' x-axis values to match on both plots. The function will fail if
#' this condition is not respected.
#'
#' @inheritParams pvalueGrob
#' @inheritParams plot_transcripts
#' @param xscale A GRanges object used to restrict the plotting region.
#'   Typically, it will be a GWAS signal or a gene known to be involved
#'   in the phenotype studied.
#' @param xexpand A numeric of length 2 representing expansion factors
#'   of x-scale limits. The first value is the expansion factor to
#'   the left while the second value is the expansion factor to the right.
#'   Values of 0 (the default) represent no expansion relative to the
#'   values in x-scale. Limits will be rounded to the nearest integer
#'   because they need to be represented as a GRanges object.
#'
#' @return NULL, invisibly. The function is called for its side-effect of
#'   plotting.
#' @export
#'
#' @examples
#' NULL
pvalue_tx_plot <- function(gwas_results, genes, transcripts, exons, cds, xscale,
			   feature = NULL,
			   xexpand = c(0, 0),
			   transcript_margins = c(0, 4.1, 0, 2.1),
			   pvalue_margins = c(5.1, 4.1, 4.1, 2.1),
			   yexpand = c(0, 0)) {

	# Checking that the left and right margins of both plots are the same number
	stopifnot(transcript_margins[2] == pvalue_margins[2] && transcript_margins[4] == pvalue_margins[4])

	# Checking that the xscale is a single range
	stopifnot(length(xscale) == 1)

	# Expanding the xscale (if applicable)
	stopifnot(length(xexpand) == 2)

	# Checking if the feature (if provided) overlaps the xscale
	if(!is.null(feature) && !overlapsAny(feature, xscale)) {
		stop("The feature parameter must overlap xscale")
	}
	
	if(!all(xexpand == 0)) {
		range_width <- GenomicRanges::width(xscale)
		GenomicRanges::start(xscale) <- round(GenomicRanges::start(xscale) - xexpand[1] * range_width)
		GenomicRanges::end(xscale) <- round(GenomicRanges::end(xscale) + xexpand[2] * range_width)
	}

	# Preparing the layout of the plot
	grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 2, heights = grid::unit(c(0.2, 0.8), "npc"))))

	# Plotting the transcripts in the top viewport
	grid::pushViewport(grid::viewport(layout.pos.row = 1))
	plot_transcripts(genes, transcripts, exons, cds, xscale = xscale, transcript_margins = transcript_margins)
	grid::upViewport()

	# Plotting the p-values in the bottom viewport
	grid::pushViewport(grid::viewport(layout.pos.row = 2))
	grid::grid.draw(pvalueGrob(gwas_results, interval = xscale,
				   feature = feature, 
				   pvalue_margins = pvalue_margins,
				   yexpand = yexpand))
	grid::upViewport()

	# Getting back up to the top viewport
	upViewport()

	invisible(NULL)
}

