# A function that generates a Manhattan plot of the results given k-mer positions, p-values, and a .fai index file
# sam_df: A sam-like data.frame of reads to which significant k-mers matched, as returned by the function add_pvalues
# fai_file: A .fai index file that is used to read chromosome sizes from
# pattern: a regular expression used to restrict the chromosomes that will be displayed on the Manhattan plot

#' kmer_manhattan
#' 
#' @param formatted_df To be completed
#' @param fai_file To be completed
#' @param pattern To be completed
#' @param signals To be completed
#'
#' @return To be completed
#'
#' @export
#' @examples
#' NULL
kmer_manhattan <- function(formatted_df, fai_file, pattern = NULL, signals = NULL) {
	fai_info <- parse_fai(fai_file, pattern)

	if(!is.null(signals) && nrow(signals)) {
		signals$first_plot <- fai_info$start[signals$first_chrom] + signals$first_pos
		signals$last_plot <- fai_info$start[signals$last_chrom] + signals$last_pos
		signals$max_plot <- fai_info$start[signals$max_chrom] + signals$max_pos
	}

	# Plotting the results
	output <- ggplot2::ggplot(formatted_df) +
		ggplot2::geom_point(mapping = ggplot2::aes(x = manhattan_rpos, y = manhattan_log10p, color = MAPQ)) +
		ggplot2::scale_x_continuous(name = "Chromosome",
					    breaks = fai_info$label_pos, 
					    labels = names(fai_info$start),
					    limits = c(0, sum(fai_info$lengths))) +
		ggplot2::ylab("-log10(p-value)") +
		ggplot2::theme_bw()

	if(!is.null(signals) && nrow(signals)) {
		output <- output + ggplot2::geom_vline(data = signals, aes(xintercept = max_plot))
		output <- output + ggplot2::geom_rect(data = signals,
						      ymin = 0, ymax = max(formatted_df$manhattan_log10p),
						      aes(xmin = first_plot, xmax = last_plot),
						      color = "red", alpha = 0.4)
	}

	return(output)
}

#' gapit_manhattan
#' 
#' @param formatted_gapit To be completed
#' @param fai_file To be completed
#' @param pattern To be completed
#' @param signals To be completed
#' @param single_panel To be completed
#' @param threshold To be completed
#'
#' @return To be completed
#'
#' @export
#' @examples
#' NULL
gapit_manhattan <- function(formatted_gapit, fai_file, pattern = NULL, signals = NULL, single_panel = TRUE, threshold = 5) {
	fai_info <- parse_fai(fai_file, pattern)

	if(!is.null(signals) && nrow(signals)) {
		signals$first_plot <- fai_info$start[signals$first_chrom] + signals$first_pos
		signals$last_plot <- fai_info$start[signals$last_chrom] + signals$last_pos
		signals$max_plot <- fai_info$start[signals$max_chrom] + signals$max_pos
	}

	# If single_panel is TRUE, then all the data is plotted along a single axis
	if(single_panel) {
		output <- ggplot2::ggplot(formatted_gapit) +
			geom_point(mapping = aes(x = manhattan_rpos, y = manhattan_log10p, color = manhattan_even)) +
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
		output <- ggplot2::ggplot(formatted_gapit) +
			geom_point(mapping = aes(x = Position, y = manhattan_log10p)) +
			facet_wrap(~manhattan_chrom, ncol = 5) +
			scale_x_continuous(name = "Position (bp)") +
			ylab("-log10(p-value)") +
			geom_hline(yintercept = threshold, linetype = 2) +
			theme_bw() +
			theme(panel.grid.minor = element_blank())

		if(!is.null(signals) && nrow(signals)) {
			output <- output + ggplot2::geom_vline(data = signals, aes(xintercept = max_plot))
			output <- output + ggplot2::geom_rect(data = signals,
							      ymin = 0, ymax = max(formatted_gapit$manhattan_log10p),
							      aes(xmin = first_pos, xmax = last_pos),
							      color = "red", alpha = 0.4)
		}
	}

	return(output)
}

