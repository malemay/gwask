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
#'
#' @return A gTree with the proper viewports and grobs set for plotting all
#'   transcripts in the provided genomic region.
#'
#' @export
#' @examples
#' NULL
transcriptsGrob <- function(genes, transcripts, exons, cds, xscale,
			    first_tx_only = FALSE) {

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
	genes <- IRanges::subsetByOverlaps(genes, xscale, ignore.strand = TRUE)

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
	main_viewport <- grid::viewport(layout = grid::grid.layout(nrow = nrows), name = "main")

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

	# Adding a rectGrob and xaxisGrob to the gTree in order to visualize the limits of the viewport
	output_gtree <- grid::addGrob(output_gtree, grid::rectGrob(gp = gpar(fill = NA)))
	output_gtree <- grid::addGrob(output_gtree, grid::xaxisGrob(vp = paste0("tx", nrows)))

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
#' @param yexpand A numeric of length 2 representing expansion factors
#'   of y-scale limits. The first value is the expansion factor to
#'   the top while the second value is the expansion factor to the bottom.
#'   Values of 0 (the default) represent no expansion relative to the
#'   values in x-scale.
#' @param merging_gap A numeric value specifying how far apart various
#'   signals in the interval GRanges object can be in order to be merged.
#'   This should be a very high value (default = 10^10) such that all ranges
#'   located on a given chromosome will be merged.
#' @param cex.points A numeric value. The expansion factor for the points.
#'
#' @return A gTree object that can be plotted with grid.draw() to display the p-values.
#'
#' @export
#' @examples
#' NULL
pvalueGrob <- function(gwas_results, interval, feature = NULL,
		       yexpand = c(0, 0), merging_gap = 10^9,
		       cex.points = 0.5) {

	if(!inherits(gwas_results, "GRanges")) {
		gwas_results <- GenomicRanges::makeGRangesFromDataFrame(gwas_results,
									keep.extra.columns = TRUE,
									ignore.strand = TRUE,
									seqnames.field = "manhattan_chrom",
									start.field = "manhattan_cpos",
									end.field = "manhattan_cpos")
	}

	# We keep only the part of the gwas_results that overlaps the interval
	gwas_results <- IRanges::subsetByOverlaps(gwas_results, interval, ignore.strand = TRUE)

	# Also adding a column with a numerical index into the signal that a marker belongs to (for coloring the points)
	if(length(gwas_results)) {
		overlaps <- IRanges::findOverlaps(gwas_results, interval)
		gwas_results$signal <- 1
		gwas_results[IRanges::from(overlaps)]$signal <- IRanges::to(overlaps)
	}

	# We also need to merge the ranges in interval so it can be used to set the x-scale for plotting
	xrange <- GenomicRanges::reduce(interval, min.gapwidth = merging_gap, ignore.strand = TRUE)
	stopifnot(length(xrange) == 1)

	# Checking if the feature (if provided) overlaps the interval
	if(!is.null(feature) && !overlapsAny(feature, xrange)) {
		warning("The feature parameter does not overlap the plotting interval; extending plotting range to include feature")
		xrange <- GenomicRanges::reduce(c(feature, xrange), min.gapwidth = merging_gap, ignore.strand = TRUE)
	}

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
	pvalue_viewport <- grid::viewport(xscale = c(GenomicRanges::start(xrange),
						     GenomicRanges::end(xrange)),
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
						       y = grid::unit(gwas_results$log10p, "native"),
						       pch = 20,
						       gp = gpar(cex = cex.points, col = gg_hue(max(gwas_results$signal))[gwas_results$signal])))

	# Adding a y-axis
	output_gtree <- grid::addGrob(output_gtree, grid::yaxisGrob())

	# Adding the axis labels
	output_gtree <- grid::addGrob(output_gtree,
				      grid::textGrob("-log10(p-value)", x = grid::unit(-3, "lines"), rot = 90))

	return(output_gtree)
}

#' Arrange a transcript grob and several p-value grobs in the same plot
#'
#' @inheritParams transcriptsGrob
#' @param pvalue_grobs A list of p-value grobs that will be arrange in
#'   viewport rows
#' @param xrange A GRanges object used to restrict the plotting region.
#'   Typically, it will be a GWAS signal or a gene known to be involved
#'   in the phenotype studied. If NULL (default), the scale is set automatically
#'   to encompass all the values provided in the pvalue_grobs.
#' @param xchrom A character of length one indicating what chromosome is being
#'   plotted. This parameter must be passed if xrange is NULL in order for
#'   the GRanges object to use for subsetting genes to be set.
#' @param tx_fraction A numeric value between 0 and 1 representing the
#'   proportion of the plot reserved for plotting the transcript(s).
#' @param margins The margins used for the plotting viewport. Default value
#'   is c(5.1, 4.1, 4.1, 2.1) (the default for \code{\link[grid]{plotViewport}})
#'
#' @return NULL, invisibly. The function is called for its side-effect of
#'   plotting.
#' @export
#'
#' @examples
#' NULL
pvalue_tx_grob <- function(pvalue_grobs, xrange = NULL, xchrom = NULL,
			   genes, transcripts, exons, cds,
			   first_tx_only = FALSE, tx_fraction = 0.1,
			   margins = c(5.1, 4.1, 4.1, 2.1)) {

	# Setting xrange and xscale for plotting
	if(is.null(xrange)) {
		if(is.null(xchrom) || !(length(xchrom) == 1)) stop("xchrom must be set to generate the GRanges objcet for subsetting genes")

		# Computing the	limits of the viewport from the limits of those in the pvalue grobs
		vp_x <- unlist(lapply(pvalue_grobs, function(x) x$vp$xscale))

		if(!is.null(vp_x)) {
			xscale <- range(vp_x)
			xrange <- GenomicRanges::GRanges(seqnames = xchrom,
							 ranges = IRanges::IRanges(start = xscale[1], end = xscale[2]))

			# Editing the x-scales of the viewports underlying pvalue_grobs so they are plotted with update coordinates
			for(i in 1:length(pvalue_grobs)) {
				if(!is.null(pvalue_grobs[[i]]$vp)) pvalue_grobs[[i]]$vp <- grid::editViewport(pvalue_grobs[[i]]$vp, xscale = xscale)
			}

		} else {	
			xrange <- NULL
		}

	} else {
		if(!is.null(xchrom)) warning("xchrom parameter ignored as xrange was set")
	}

	# Preparing the layout of the plot
	# There is one row for the transcripts, and one row per pvalue grob
	nrows <- length(pvalue_grobs) + 1
	main_viewport <- grid::viewport(layout = grid::grid.layout(nrow = nrows, heights = unit(c(tx_fraction,  rep(1, nrows - 1)),
												c("npc", rep("null", nrows - 1)))))

	# Creating all the row viewports
	sub_viewports <- lapply(1:nrows, function(x) grid::vpStack(grid::viewport(layout.pos.row = x,
										  name = paste0("row_", x)),
								   grid::plotViewport(margins = if(x != 1) margins else c(0, margins[2], 0, margins[4]),
										      name = "plotting_region")))

	# Creating the output gTree object with the main viewport and childrenvp set
	output_gtree <- grid::gTree(vp = main_viewport,
				    childrenvp = do.call("vpList", sub_viewports),
				    cl = "pvalue_tx_grob")

	# Creating a gTree that will contain all the transcripts to be plotted
	# This is only done if xrange is not NULL
	if(!is.null(xrange)) {
		# Plotting the transcripts in the top viewport
		tx_gtree <- grid::gTree(children = grid::gList(transcriptsGrob(genes, transcripts, exons, cds,
									       first_tx_only = first_tx_only,
									       xscale = xrange)),
					vp = grid::vpPath("row_1", "plotting_region"))

		output_gtree <- grid::addGrob(output_gtree, tx_gtree)
	}

	# Plotting the p-values in the bottom viewport

	for(i in 1:length(pvalue_grobs)) {
		output_gtree <- grid::addGrob(output_gtree,
					      grid::gTree(children = gList(pvalue_grobs[[i]]),
							  vp = vpPath(paste0("row_", as.character(i + 1)), "plotting_region")))
	}


	return(output_gtree)
}

#' Generating a vector of colors comparable to those used by ggplot2
#'
#' @param n An integer of length one. The number of colors to generate
#'
#' @return A vector of colors to be used for plotting
#'
#' @examples
#' NULL
gg_hue <- function(n) {
	hues = seq(15, 375, length = n + 1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}

