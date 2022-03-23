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
#'
#' @return A ggplot object describing a Manhattan plot with the specified parameters.
#'
#' @export
#' @examples
#' NULL
manhattan_plot <- function(formatted_data, gwas_type, fai_file, pattern = NULL,
			   signals = NULL, extent = FALSE, single_panel = TRUE,
			   threshold = 5, min_log10p = 0) {

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
									  color = ifelse(gwas_type == "gapit",
											 "manhattan_even",
											 "MAPQ"))) +
				ggplot2::scale_x_continuous(name = "Chromosome",
							    breaks = fai_info$label_pos,
							    labels = names(fai_info$start)) +
				ggplot2::ylab("-log10(p-value)") +
				ggplot2::geom_hline(yintercept = threshold) +
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
			ggplot2::geom_hline(yintercept = threshold, linetype = 2) +
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

	return(output)
}

