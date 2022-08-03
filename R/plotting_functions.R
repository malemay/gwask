#' A function that returns a graphical object (grob) representing a manhattan plot
#' 
#' @param gwas_results A GRanges object formatted with either \code{\link{format_kmer_gwas}}
#'   or \code{\link{format_gapit_gwas}}. These functions
#'   return GRanges containing standard columns and seqinfo in addition to GWAS type-specific
#'   columns, therefore allowing a common interface for all GWAS types.
#' @param threshold Numeric. A single value indicating where the -log10(p-value)
#'   threshold line should be plotted.
#' @param min_log10p Numeric. The minimum -log10(p-value) to display on the plot
#'   can be used to filter out some non-significant markers for plotting speedup.
#' @param signals A GenomicRanges object containing reference signals to be added
#'   to the plot. If NULL (the default), then no signals are plotted.
#' @param numeric_chrom Logical. Whether chromosome names should be stripped from
#'   their alphabetical component and converted to numeric values (default: FALSE).
#' @param yexpand A numeric of length two. A fractional value relative to the
#'   range of the y-axis used to expand the y-scale on either side. The first
#'   value is the expansion factor at the bottom of the scale, and the second
#'   value is the expansion factor at the top of the scale.
#' @param xexpand Same as yexpand but for the x-scale.
#' @param cex.points A numeric value. The expansion factor for the points.
#' @param cex.lab A numeric value. The expansion factor for the axis labels.
#' @param margins The margins used for the plotting viewport. Default value
#'   is c(5.1, 4.1, 4.1, 2.1) (the default for \code{\link[grid]{plotViewport}})
#'
#' @return A gTree describing a Manhattan plot with the specified parameters.
#'
#' @export
#' @examples
#' NULL
manhattanGrob <- function(gwas_results, threshold = NULL, min_log10p = 0,
			  signals = NULL, numeric_chrom = FALSE,
			  yexpand = c(0.03, 0.08), xexpand = c(0.03, 0.03),
			  cex.points = 0.5, cex.lab = 1,
			  margins = c(5.1, 4.1, 4.1, 2.1)) {

	# Checking that gwas_results has the right format
	if(!"log10p" %in% names(GenomicRanges::mcols(gwas_results))) stop("The required log10p column is missing from gwas_results.")

	stopifnot(length(yexpand) == 2)
	stopifnot(length(xexpand) == 2)

	# Using the seqinfo of the gwas_results input to compute various information
	chrom_lengths <- GenomeInfoDb::seqlengths(gwas_results)

	# A vector of the starting positions of chromosomes on a linear scale
	chrom_start <- cumsum(chrom_lengths)
	chrom_start <- c(0, chrom_start[-(length(chrom_start))]) + 1

	# A vector of the positions where the x-axis labels should be drawn on a linear scale
	label_pos <- chrom_start + 0.5 * chrom_lengths

	# All vectors have the names of the sequences to allow lookup of the info by name
	names(chrom_start) <- names(chrom_lengths) <- names(label_pos) <- GenomeInfoDb::seqlevels(gwas_results)

	# Computing the absolute position on the genome as the median between start and end pos
	gwas_results$manhattan_rpos <- (GenomicRanges::start(gwas_results) + GenomicRanges::end(gwas_results)) / 2 +
		chrom_start[as.character(GenomicRanges::seqnames(gwas_results))]

	# Computing whether the chromosome chromsome is even or odd (for plotting purposes)
	even_chromosomes <- GenomeInfoDb::seqlevels(gwas_results)[c(FALSE, TRUE)]
	gwas_results$manhattan_even <- as.character(GenomicRanges::seqnames(gwas_results)) %in% even_chromosomes

	# Removing the associations that are less significant than min_log10p
	if(min_log10p > 0) gwas_results <- gwas_results[gwas_results$log10p >= min_log10p]

	# If signals are provided, we need to also transform their coordinates to linear positions
	if(!is.null(signals)) {
		signals$manhattan_rpos <- (GenomicRanges::start(signals) + GenomicRanges::end(signals)) / 2 +
			chrom_start[as.character(GenomicRanges::seqnames(signals))]
		# Also rescaling the p-values to npc coordinates between 0.2 and 0.8 for plotting
		signals$pvalue_rescaled <- signals$log_pvalue / max(signals$log_pvalue) * 0.6 + 0.2
	}

	# Setting the plotting scales
	if(length(gwas_results)) {
		yscale <- c(0, max(gwas_results$log10p, threshold))

		if(!all(yexpand == 0)) {
			yrange <- diff(yscale)
			yscale[1] <- yscale[1] - yexpand[1] * yrange
			yscale[2] <- yscale[2] + yexpand[2] * yrange
		}
	} else {
		yscale <- c(0, threshold + 1)
	}

	xscale <- c(0, sum(chrom_lengths))

	if(!all(xexpand == 0)) {
		xrange <- diff(xscale)
		xscale[1] <- xscale[1] - xexpand[1] * xrange
		xscale[2] <- xscale[2] + xexpand[2] * xrange
	}

	# Initializing the output gTree object
	output_gtree <- grid::gTree(vp = grid::plotViewport(margins = margins,
							    xscale = xscale,
							    yscale = yscale,
							    name = "manhattan_vp"),
				    name = "manhattan")

	# Adding a box around the viewport
	output_gtree <- grid::addGrob(output_gtree, rectGrob(name = "manahttan_box"))

	# Adding the points
	if(length(gwas_results)) {
		output_gtree <- grid::addGrob(output_gtree,
					      grid::pointsGrob(x = gwas_results$manhattan_rpos,
							       y = gwas_results$log10p,
							       pch = 20,
							       gp = grid::gpar(cex = cex.points, col = ifelse(gwas_results$manhattan_even, "red", "blue")),
							       name = "manhattan_points"))
	}

	# Adding the x-axis
	output_gtree <- grid::addGrob(output_gtree,
				      grid::xaxisGrob(at = label_pos,
						      label = if(numeric_chrom) as.character(as.numeric(gsub("[a-zA-Z]", "", names(label_pos)))) else names(label_pos),
						      gp = grid::gpar(cex = cex.lab),
						      name = "manhattan_xaxis"))

	# Adding the x-axis label
	output_gtree <- grid::addGrob(output_gtree, grid::textGrob("Chromosome", y = grid::unit(-3, "lines"), name = "manhattan_xlabel"))

	# Adding the y-axis
	output_gtree <- grid::addGrob(output_gtree, grid::yaxisGrob(gp = grid::gpar(cex = cex.lab), name = "manhattan_yaxis"))

	# Adding the x-axis label
	output_gtree <- grid::addGrob(output_gtree, grid::textGrob("-log10(p-value)", x = grid::unit(-3, "lines"), rot = 90, name = "manhattan_ylabel"))

	# Also adding the threshold if not NULL
	if(!is.null(threshold)) {
		output_gtree <- grid::addGrob(output_gtree,
					      grid::linesGrob(x = c(0, 1),
							      y = grid::unit(threshold, "native"),
							      gp = grid::gpar(lty = "13"),
							      name = "manhattan_threshold"))
	}

	# Adding the signals if not NULL
	if(!is.null(signals)) {
		for(i in 1:length(signals)) {
			output_gtree <- grid::addGrob(output_gtree,
						      grid::linesGrob(x = grid::unit(signals[i]$manhattan_rpos, "native"),
								      y = c(0, 1),
								      gp = grid::gpar(lty = "16"),
								      name = paste0("manhattan_sigline_", i)))

			output_gtree <- grid::addGrob(output_gtree,
						      grid::textGrob(label = signals[i]$locus,
								     x = grid::unit(signals[i]$manhattan_rpos, "native"),
								     y = signals[i]$pvalue_rescaled, hjust = 0.1,
								     name = paste0("manhattan_siglabel_", i)))
		}
	}

	return(output_gtree)
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
#' \code{\link{transcriptsGrob}} that will take care of organizing
#' the viewports and coordinate scales for plotting all the transcripts
#' of a given gene model.
#'
#' @param tx_gene A GRanges object representing the coordinates of
#'   of a single gene to which the transcript belongs.
#' @param tx_exons A GRanges object representing the exons belonging
#'   to a single transcript.
#' @param tx_name A GRanges object representing the coding sequences
#'   belonging to a single transcript.
#' @param output_vp A viewport that the transcript should be plotted in.
#' @param name The name of the returned gTree object.
#'
#' @return A gTree object representing a single transcript.
#'
#' @export
#' @examples
#' NULL
transcriptGrob <- function(tx_gene, tx_exons, tx_cds, output_vp = NULL, name = NULL) {

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

	return(gTree(name = name,
		     cl = "txgrob",
		     children = gList(tx_line, tx_exons, tx_cds),
		     vp = output_vp))
}

#' Generate a grob of all the possible transcripts for genes in a genomic region using grid
#'
#' This function calls transcriptGrob on all individual transcripts
#' of the genes located in a given genomic region and organizes each separately
#' in a row of a viewport grid. The object returned is a gTree that can be plotted
#' with grid.draw().
#'
#' This function requires as input the full set of genes, transcripts,
#' exons and coding sequences for the reference genome being used. It 
#' will check that transcripts exist for all genes provided. This
#' behaviour may be modified in the future depending on use cases.
#'
#' It is also an error in the current implementation to try and plot
#' over an interval that does not overlap any genes.
#'
#' The function itself does not add any grobs by itself
#' but instead arranges the viewport layout and calls the
#' \code{\link{transcriptGrob}} to add the transcripts one by one.
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
#' @param first_tx_only A logical. Whether to plot only the first transcript
#'   of a given gene (default: FALSE).
#' @param xscale A GRanges object used to set the limits of the x-axis
#'   and to subset the genes object such that only the ones located in the 
#'   interval are plotted.
#' @param transcript_margins A numeric of length 4 or that can be recycled to that
#'   length. The margins are used by the function \code{\link[grid]{plotViewport}}
#'   to set the margins around the plotting region. By default it is set
#'   to c(0, 4.1, 0, 2.1).
#'
#' @return A gTree with the proper viewports and grobs set for plotting all
#'   transcripts in the provided genomic region.
#'
#' @export
#' @examples
#' NULL
transcriptsGrob <- function(genes, transcripts, exons, cds, xscale,
			    first_tx_only = FALSE,
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

	if(!length(genes)) {
		warning("No genes in the selected interval")
		return(gTree())
	}

	# Extracting the transcripts related to those genes
	transcripts <- transcripts[names(genes)]

	# The number of rows for the plotting viewports depends on the gene with the most transcripts
	nrows <- max(lengths(transcripts))
	if(first_tx_only) nrows <- 1

	# A viewport layout with as many rows as there are transcripts for the gene with the highest number of transcripts
	main_viewport <- grid::plotViewport(margins = transcript_margins,
					    layout = grid::grid.layout(nrow = nrows),
					    name = "main")

	# Creating all the row viewports
	sub_viewports <- lapply(1:nrows, function(x) grid::viewport(layout.pos.row = x,
								    xscale = c(GenomicRanges::start(xscale),
									       GenomicRanges::end(xscale)),
								    name = paste0("tx", x)))

	# Creating a gTree that will contain all the transcripts to be plotted
	output_gtree <- gTree(vp = main_viewport,
			      childrenvp = do.call("vpList", sub_viewports),
			      cl = "txxgrob")

	# Iterating over the genes
	for(gene_index in 1:length(genes)) {

		# Getting the names of the transcripts of that gene
		tnames <- transcripts[[ names(genes)[gene_index] ]]$tx_name
		if(first_tx_only) tnames <- tnames[1]

		# Looping over all transcripts to plot them in separate viewports
		for(i in 1:length(tnames)) {
			tx_name <- tnames[i]

			# Drawing a transcriptGrob for each transcript
			output_gtree <- grid::addGrob(output_gtree,
						      transcriptGrob(genes[gene_index],
								     exons[[tx_name]],
								     cds[[tx_name]],
								     name = paste0("gene", gene_index, "_tx", i),
								     output_vp = vpPath(paste0("tx", i))))
		}
	}

	return(output_gtree)
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
		# If the y-scale interval is only one marker, we set its minimum to 0
		if(yscale[1] == yscale[2]) yscale[1] <- 0
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
#' @inheritParams transcriptsGrob
#' @param xscale A GRanges object used to restrict the plotting region.
#'   Typically, it will be a GWAS signal or a gene known to be involved
#'   in the phenotype studied.
#' @param xexpand A numeric of length 2 representing expansion factors
#'   of x-scale limits. The first value is the expansion factor to
#'   the left while the second value is the expansion factor to the right.
#'   Values of 0 (the default) represent no expansion relative to the
#'   values in x-scale. Limits will be rounded to the nearest integer
#'   because they need to be represented as a GRanges object.
#' @param pvalue_fraction A numeric value between 0 and 1. The proportion
#'   of the plotting viewport reserved for the p-value grob. Default: 0.8.
#'
#' @return NULL, invisibly. The function is called for its side-effect of
#'   plotting.
#' @export
#'
#' @examples
#' NULL
pvalue_tx_grob <- function(gwas_results, genes, transcripts, exons, cds, xscale,
			   first_tx_only = FALSE,
			   feature = NULL,
			   xexpand = c(0, 0),
			   transcript_margins = c(0, 4.1, 0, 2.1),
			   pvalue_margins = c(5.1, 4.1, 4.1, 2.1),
			   yexpand = c(0, 0),
			   pvalue_fraction = 0.8) {

	# Checking that the left and right margins of both plots are the same number
	stopifnot(transcript_margins[2] == pvalue_margins[2] && transcript_margins[4] == pvalue_margins[4])

	# Checking that the xscale is a single range
	stopifnot(length(xscale) == 1)

	# Expanding the xscale (if applicable)
	stopifnot(length(xexpand) == 2)

	# Checking if the feature (if provided) overlaps the xscale
	if(!is.null(feature) && !overlapsAny(feature, xscale)) {
		warning("The feature parameter does not overlap xscale; setting feature to NULL")
		feature = NULL
	}
	
	if(!all(xexpand == 0)) {
		range_width <- GenomicRanges::width(xscale)
		GenomicRanges::start(xscale) <- round(GenomicRanges::start(xscale) - xexpand[1] * range_width)
		GenomicRanges::end(xscale) <- round(GenomicRanges::end(xscale) + xexpand[2] * range_width)
	}

	# Preparing the layout of the plot
	main_viewport <- grid::viewport(layout = grid::grid.layout(nrow = 2, heights = grid::unit(c(1, pvalue_fraction), c("null", "npc"))))

	output_gtree <- grid::gTree(vp = main_viewport,
				    childrenvp = vpList(grid::viewport(layout.pos.row = 1, name = "tx_vp"),
							grid::viewport(layout.pos.row = 2, name = "pvalue_vp")))

	# Plotting the transcripts in the top viewport
	tx_gtree <- grid::gTree(children = gList(transcriptsGrob(genes, transcripts, exons, cds,
								 first_tx_only = first_tx_only,
								 xscale = xscale,
								 transcript_margins = transcript_margins)),
				vp = "tx_vp")

	# Plotting the p-values in the bottom viewport
	pvalue_gtree <- grid::gTree(children = gList(pvalueGrob(gwas_results, interval = xscale,
								feature = feature, 
								pvalue_margins = pvalue_margins,
								yexpand = yexpand)),
				    vp = "pvalue_vp")

	output_gtree <- addGrob(output_gtree, tx_gtree)
	output_gtree <- addGrob(output_gtree, pvalue_gtree)

	return(output_gtree)
}

