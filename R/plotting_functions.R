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


