#' manhattan_plot
#' 
#' @param formatted_data A data.frame formatted with either \code{\link{format_kmer_gwas}}
#'   or \code{\link{format_gapit_gwas}} depending on the GWAS type. These functions
#'   return data.frames containing standard columns in addition to GWAS type-specific
#'   columns, therefore allowing a common interface for all GWAS types.
#' @param gwas_type A character of length one. Only two options are supported at
#'   this time, either \code{"kmer"} for k-mer GWAS as implemented by Voichek
#'   and Weigel (2020) or \code{"gapit"} for the output of GAPIT.
#' @param fai_file A character of length one. The name of the .fai index file
#'   to extract the reference genome data from.
#' @param pattern A character of length one. A regular expression to use
#'   for selecting what reference sequences to keep. If \code{NULL}
#'   (the default), all sequences are used. The function will force all
#'   references sequences kept to be plotted even if no data point or
#'   signal is observed on it.
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
manhattan_plot <- function(formatted_data, gwas_type, fai_file, pattern = NULL,
			   signals = NULL, extent = FALSE, single_panel = TRUE,
			   threshold = NULL, min_log10p = 0, plot_mapq = FALSE) {

	# Sanity checks -------------------------------------

	# Checking that formatted_data has the right format
	if(!all(c("manhattan_chrom", "manhattan_rpos", "manhattan_cpos", "manhattan_log10p", "manhattan_even") %in% names(formatted_data))) {
		stop("The formatted_data data.frame does not satisfy format requirements.")
	}

	# Checking that gwas_type is supported
	if(!gwas_type %in% c("kmer", "gapit")) {
		stop("gwas_type ", gwas_type, " is not supported.")
	}

	# Checking that the MAPQ column is suppled if type == kmer
	if(gwas_type == "kmer" && ! "MAPQ" %in% names(formatted_data)) {
		stop("MAPQ column must be present in formatted_data for gwas_type = \"kmer\"")
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
	if(plot_mapq && ! "MAPQ" %in% colnames(formatted_data)) {
		warning("plot_mapq = TRUE but \"MAPQ\" is not among the columns. Setting plot_mapq to FALSE")
		plot_mapq <- FALSE
	}

	# End of sanity checks ------------------------------

	# Parsing the .fai index file to extract information about the reference sequences
	fai_info <- parse_fai(fai_file, pattern)

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
		formatted_data <- formatted_data[formatted_data$manhattan_log10p >= min_log10p, ]
	}

	# If single_panel is TRUE, then all the data is plotted along a single axis
	if(single_panel) {
		output <- ggplot2::ggplot(formatted_data) +
			ggplot2::geom_point(mapping = ggplot2::aes_string(x = "manhattan_rpos",
									  y = "manhattan_log10p",
									  color = ifelse(plot_mapq, "MAPQ", "manhattan_even"))) +
				ggplot2::scale_x_continuous(name = "Chromosome",
							    breaks = fai_info$label_pos,
							    labels = names(fai_info$start),
							    limits = c(0, sum(fai_info$lengths))) +
				ggplot2::ylab("-log10(p-value)") +
				ggplot2::theme_bw() +
				ggplot2::theme(panel.grid.minor = ggplot2::element_blank())


	} else {
	# If single_panel is not TRUE, there is one panel per chromosome
		output <- ggplot2::ggplot(formatted_data) +
			ggplot2::geom_point(mapping = ggplot2::aes_string(x = "Position",
									  y = "manhattan_log10p")) +
			ggplot2::facet_wrap(~manhattan_chrom, ncol = 5) +
			ggplot2::scale_x_continuous(name = "Position (bp)") +
			ggplot2::ylab("-log10(p-value)") +
			ggplot2::theme_bw() +
			ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
	}

	if(!is.null(signals) && nrow(signals)) {
		output <- output + ggplot2::geom_vline(data = signals, aes(xintercept = max_plot), linetype = 2)

		if(extent) {
			output <- output + ggplot2::geom_rect(data = signals,
							      ymin = 0, ymax = max(formatted_data$manhattan_log10p),
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

#' Plot a transcript using grid functions
#'
#' This function plots a transcript model using grid functions.
#' 
#' The function assumes that that plotting occurs in a viewport
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
#' @return NULL, invisibly. The function is invoked for its side
#'   effect of plotting the transcript.
#'
#' @export
#' @examples
#' NULL
plot_transcript <- function(tx_gene, tx_exons, tx_cds) {

	# Drawing the line that goes from start to end of the gene
	grid::grid.lines(x = grid::unit(c(GenomicRanges::start(tx_gene), GenomicRanges::end(tx_gene)), "native"),
			 y = grid::unit(0.5, "npc"))

	# Drawing the rectangles that correspond to exons
	grid::grid.rect(x = grid::unit(GenomicRanges::start(tx_exons), "native"),
			y = grid::unit(0.5, "npc"),
			width = grid::unit(GenomicRanges::width(tx_exons), "native"),
			height = grid::unit(0.8, "npc"),
			hjust = 0,
			gp = grid::gpar(fill = "white"))

	# Doing the same thing for the coding sequences by using a shaded color
	grid::grid.rect(x = grid::unit(GenomicRanges::start(tx_cds), "native"),
			y = grid::unit(0.5, "npc"),
			width = grid::unit(GenomicRanges::width(tx_cds), "native"),
			height = grid::unit(0.8, "npc"),
			hjust = 0,
			gp = grid::gpar(fill = ifelse(as.character(GenomicRanges::strand(tx_cds)) == "+", "skyblue", "orange")))

	return(invisible(NULL))
}

#' Plot all the possible transcripts for a gene using grid functions
#'
#' This function calls plot_transcript on all individual transcripts
#' of a gene and plots each separately in a row of a viewport grid.
#'
#' At the moment, this function requires as input the full set of
#' genes, transcripts, exons and coding sequences for the reference
#' genome being used. It will check that transcripts exist for all
#' genes provided even if it is meant to plot a single gene. This
#' behaviour may be modified in the future depending on use cases.
#'
#' The function itself does not plot any objects but instead arranges
#' the viewport layout and moving between viewports, and calls the
#' \code{\link{plot_transcript}} to plot the transcripts one by one.
#' The latter is the function that should be modified for changing
#' the way that plotted transcripts look.
#'
#' @param gene_name A character. The name of the gene whose transcripts
#'   are to be plotted.
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
plot_transcripts <- function(gene_name, genes, transcripts, exons, cds,
			     transcript_margins = c(0, 4.1, 0, 2.1)) {

	# Checking that the inputs are of the right type
	stopifnot(inherits(genes, "GRanges"))
	stopifnot(inherits(transcripts, "GRangesList"))
	stopifnot(inherits(exons, "GRangesList"))
	stopifnot(inherits(cds, "GRangesList"))

	# Doing a few more sanity checks on the inputs
	stopifnot(gene_name %in% names(genes))
	stopifnot(all(names(genes) %in% names(transcripts)))
	stopifnot(sum(lengths(transcripts)) == length(exons))
	stopifnot(sum(lengths(transcripts)) == length(cds))

	# Extracting the gene itself from the genes object
	gene <- genes[gene_name]

	# Getting the names of the transcripts of that gene
	tnames <- transcripts[[gene_name]]$tx_name

	# Pushing a viewport layout with as many rows as there are transcripts for that gene
	grid::pushViewport(grid::plotViewport(margins = transcript_margins, layout = grid::grid.layout(nrow = length(tnames))))

	# Looping over all transcripts to plot them in separate viewports
	for(i in 1:length(tnames)) {
		tx_name <- tnames[i]
		grid::pushViewport(grid::plotViewport(layout.pos.row = i,
						      xscale = c(GenomicRanges::start(gene),
								 GenomicRanges::end(gene))))

		# Calling the plot_transcript function for each transcript
		plot_transcript(gene, exons[[tx_name]], cds[[tx_name]])

		# Moving back to the parent viewport
		grid::upViewport()
	}

	# Going back to the viewport that we started from
	grid::upViewport()

	return(invisible(NULL))
}

#' Plot the p-values of a GWAS analysis at a locus using grid functions
#'
#' This function takes a data.frame of formatted GWAS results and plots
#' the p-values in a given interval.
#'
#' @param gwas_results A GRanges object or data.frame of GWAS results.
#'   If a data.frame, must include the columns manhattan_chrom and
#'   manhattan_cpos for coercion to a GRanges. The metadata column
#'   manhattan_log10p must be present as it will be used to plot the
#'   p-values on the y-axis.
#' @param interval A GRanges object used to subset the plotting to a given
#'   region.
#' @param pvalue_margins A numeric of length 4 or that can be recycled to that
#'   length. The margins used by the function \code{\link[grid]{plotViewport}}
#'   to set the margins around the plotting region. By default it is set
#'   to c(5.1, 4.1, 4.1, 2.1).
#'
#' @return NULL, invisibly. This function is called for its plotting
#'   side-effects.
#'
#' @export
#' @examples
#' NULL
pvalue_plot <- function(gwas_results, interval, pvalue_margins  = c(5.1, 4.1 , 4.1, 2.1)) {
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

	# It is an error to try and plot empty data
	if(!length(gwas_results)) stop("No GWAS data in interval")

	# Otherwise we can go on with the plotting by creating the viewport with appropriate scales
	grid::pushViewport(grid::plotViewport(margins = pvalue_margins,
					      xscale = c(GenomicRanges::start(interval),
							 GenomicRanges::end(interval)),
					      yscale = c(min(gwas_results$manhattan_log10p),
							 max(gwas_results$manhattan_log10p))))

	# We plot a box around the viewport
	grid.rect()

	# Then we can plot the data points
	grid::grid.points(x = grid::unit(GenomicRanges::start(gwas_results), "native"),
			  y = grid::unit(gwas_results$manhattan_log10p, "native"))

	# Adding x- and y-axis
	grid.xaxis()
	grid.yaxis()

	# Adding the axis labels
	grid.text("-log10(p-value)", x = grid::unit(-3, "lines"), rot = 90)
	grid.text("Position along reference (bp)", y = grid::unit(-3, "lines"))

	# Going back up one viewport to where we were before
	grid::upViewport()
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
#' @inheritParams pvalue_plot
#' @inheritParams plot_transcripts
#'
#' @return NULL, invisibly. The function is called for its side-effect of
#'   plotting.
#' @export
#'
#' @examples
#' NULL
pvalue_tx_plot <- function(gwas_results, gene_name, genes, transcripts, exons, cds,
			   transcript_margins = c(0, 4.1, 0, 2.1),
			   pvalue_margins = c(5.1, 4.1, 4.1, 2.1)) {

	# Checking that the left and right margins of both plots are the same number
	stopifnot(transcript_margins[2] == pvalue_margins[2] && transcript_margins[4] == pvalue_margins[4])

	gene <- genes[gene_name]

	# Preparing the layout of the plot
	grid::grid.newpage()
	grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 2, heights = grid::unit(c(0.2, 0.8), "npc"))))

	# Plotting the p-values in the top viewport ; for this we need to extract the limits of the gene
	grid::pushViewport(grid::viewport(layout.pos.row = 2))
	pvalue_plot(gwas_results, gene, pvalue_margins = pvalue_margins)
	grid::upViewport()

	# Plotting the transcripts in the bottom viewport
	grid::pushViewport(grid::viewport(layout.pos.row = 1))
	plot_transcripts(gene_name, genes, transcripts, exons, cds, transcript_margins = transcript_margins)
	grid::upViewport()

	invisible(NULL)
}

