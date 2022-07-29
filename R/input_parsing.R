#' Format k-mer GWAS results for downstream analyses
#'
#' This function takes a sorted and indexed .bam file of reads that matched 
#' significant k-mers as obtained from the katcher program and annotated
#' with p-values by the add_pvalues program. It returns a GRanges object
#' formatted for downstream analyses, including manhattan plotting.
#'
#' @param input_bam A character. The path to an indexed .bam file that
#'   contains KM and PV tags indicating the matching k-mers and their respective
#'   p-values
#' @param ref_fasta A character of length one. A fasta file that
#'   will be queried with the \code{\link[Rsamtools]{scanFaIndex}} function to
#'   get information about the reference genome and set the seqinfo fields
#'   of the output GRanges object.
#' @param pattern A character of length one. A regular expression
#'   used to match seqlevels of the chromsomes or scaffolds in the
#'   fasta index. Matching sequences will be kept. If NULL
#'   (the default), all sequences are kept. In the current implementation,
#'   an error will be thrown if changing the seqlevels results in matches
#'   being pruned. NO LONGER TRUE
#' @param min_mapq A numeric of length one. The minimum mapping quality
#'   of the alignment on which the k-mer was found for it to be kept.
#'   The default value (0) implies no filtering.
#'
#' @return A GRanges object that contains the information on matching
#'   k-mers and can be used for downstream analyses including the
#'   function \code{\link{manhattan_plot}}. Must contain a field
#'   called log10_p that stores the -log10(p-value) of the k-mer.
#'
#' @export
#' @examples
#' NULL
format_kmer_gwas <- function(input_bam, ref_fasta, pattern = NULL, min_mapq = 0) {

	# Reading the .fai index info
	stopifnot(file.exists(paste0(ref_fasta, ".fai")))
	fai_info <- Rsamtools::scanFaIndex(ref_fasta)

	# Reading the .bam file
	bam_records <- Rsamtools::scanBam(input_bam,
					  param = Rsamtools::ScanBamParam(what = Rsamtools::scanBamWhat(),
									  tag = c("KM", "PV")))[[1]]

	# Generating a data.frame from the records and formatting the columns
	sam_df <- as.data.frame(bam_records[1:11])

	# We need to handle the case if there were no reads in the .bam file
	if(!nrow(sam_df)) {
		output <- GRanges()

		output$qname <- character()
		output$flag <- integer()
		output$pos <- integer()
		output$qwidth <- integer()
		output$mapq <- integer()
		output$cigar <- character()
		output$mrnm <-  character()
		output$mpos <- integer()
		output$isize <- integer()
		output$pvalue <- numeric()
		output$kmer_canon <- character()
		output$kmer_pos_count <- integer()
		output$log10p <- numeric()

		return(output)
	}

	sam_df$rname <- as.character(sam_df$rname)
	sam_df$strand <- as.character(sam_df$strand)
	sam_df$mrnm <- as.character(sam_df$mrnm)

	# Adding the columns that require particular formatting
	sam_df$seq <- as.character(bam_records$seq)
	sam_df$qual <- as.character(bam_records$qual)
	sam_df$KM <- as.character(bam_records$tag$KM)
	sam_df$PV <- as.character(bam_records$tag$PV)

	# Removing the alignments for which the reference chromosome or position is unknown
	sam_df <- sam_df[!is.na(sam_df$rname) & !is.na(sam_df$pos),]

	# The p-value is the highest (most significant) value of any k-mer found in the read
	# We extract both the minimum p-value and the k-mer associated with it
	min_pv <- sapply(strsplit(sam_df$PV, ","), function(x) which.min(as.numeric(x)))
	sam_df$kmer <- mapply(function(x, y) x[y], strsplit(sam_df$KM, ","), min_pv)
	sam_df$pvalue <- mapply(function(x, y) as.numeric(x[y]), strsplit(sam_df$PV, ","), min_pv)

	# Getting the position of the k-mer with the most significant value within the read
	sam_df$kpos <- sam_df$pos + mapply(function(x, y) regexpr(x, y, fixed = TRUE), sam_df$kmer, sam_df$seq)

	# We add a column for the reverse complement of the k-mer sequence, and another for the canonized version of the k-mer
	sam_df$kmer_rev   <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(sam_df$kmer)))
	sam_df$kmer_canon <- ifelse(sam_df$kmer < sam_df$kmer_rev, sam_df$kmer, sam_df$kmer_rev)

	# We remove reads with k-mers that match at the same position, but first we count their number of occurences
	sam_df$kmer_pos_label <- paste0(sam_df$rname, "_", sam_df$kpos, "_", sam_df$kmer_canon)
	kmer_pos_table <- table(sam_df$kmer_pos_label)
	sam_df$kmer_pos_count <- as.integer(kmer_pos_table[sam_df$kmer_pos_label])

	sam_df <- sam_df[!duplicated(sam_df[, c("rname", "kmer_canon", "kpos")]), ]

	# Calculating the end of the ranges from the length of the kmers
	sam_df$end <- sam_df$kpos + nchar(sam_df$kmer) - 1

	# Filtering based on the MAPQ
	sam_df <- sam_df[sam_df$mapq >= min_mapq, ]

	# Coercing the input data.frame to a GRanges object
	output <- GenomicRanges::makeGRangesFromDataFrame(sam_df,
							  start.field = "kpos",
							  seqnames.field = "rname",
							  end.field = "end",
							  keep.extra.columns = TRUE)

	# Setting the seqlevels and seqlengths of the output GRanges from the .fai index
	GenomeInfoDb::seqlevels(output) <- GenomeInfoDb::seqlevels(fai_info)
	GenomeInfoDb::seqlengths(output) <- GenomeInfoDb::seqlengths(fai_info)

	# Pruning the levels according to the pattern
	if(!is.null(pattern)) GenomeInfoDb::seqlevels(output, pruning.mode = "coarse") <- grep(pattern, GenomeInfoDb::seqlevels(output), value = TRUE)

	# We keep only the position at which a k-mer is most often observed
	output <- output[order(output$kmer_canon, output$kmer_pos_count, decreasing = TRUE)]
	output <- output[!duplicated(output$kmer_canon)]
	output <- sort(output)

	# Removing unnecessary memory-hungry columns
	output$seq <- NULL
	output$qual <- NULL
	output$KM <- NULL
	output$PV <- NULL
	output$kmer <- NULL
	output$kmer_rev <- NULL
	output$kmer_pos_label <- NULL

	# Adding a column for the -log10(p-value)
	output$log10p <- -log10(output$pvalue)

	return(output)
}

#' Format GAPIT GWAS results for generating Manhattan plots
#'
#' This function takes a CSV file of GAPIT GWAS results and formats
#' them as a GRanges object for processing with other functions, including
#' \code{\link{manhattan_plot}}.
#'
#' Only relevant fields (Chromosome, Position, P.value, FDR_Adjusted_P-values)
#' are kept to increase reading efficiency.
#'
#' @param filename Character of length one. The name of the .csv file
#'   containing the results of the GAPIT analysis.
#' @param ref_fasta A character of length one. A fasta file that
#'   will be queried with the \code{\link[Rsamtools]{scanFaIndex}} function to
#'   query information about the reference genome and set the seqinfo fields
#'   of the output GRanges object.
#' @param chromosomes a character vector of the same length
#'   as the number of unique values in the Chromosome column
#'   of the GAPIT output. The integer chromosome values in the
#'   results returned by GAPIT will be used to index into that
#'   character vector to translate integer values into chromosome
#'   names.
#' @param pattern A character of length one. A regular expression
#'   used to match seqlevels of the chromsomes or scaffolds in the
#'   fasta index. Matching sequences will be kept. If NULL
#'   (the default), all sequences are kept. In the current implementation,
#'   an error will be thrown if changing the seqlevels results in matches
#'   being pruned.
#'
#' @return A GRanges object that contains the information on associated
#'   markers and can be used for downstream analyses including the
#'   function \code{\link{manhattan_plot}}. Must contain a field
#'   called log10_p that stores the -log10(p-value) of the marker.
#'
#' @export
#' @examples
#' NULL
format_gapit_gwas <- function(filename, ref_fasta, chromosomes, pattern = NULL) {

	# Reading the file with the GAPIT results
	gapit_df <- read.table(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE,
			       colClasses = c("NULL", "integer", "numeric", "numeric", "NULL", "NULL", "NULL", "NULL", "numeric", "NULL"))

	# Reading the .fai index info
	stopifnot(file.exists(paste0(ref_fasta, ".fai")))
	fai_info <- Rsamtools::scanFaIndex(ref_fasta)

	# Checking that the number of chromosomes in gapit_df matches
	# the number of chromosomes in the chromosomes character vector
	n_chrom <- length(unique(gapit_df$Chromosome))
	stopifnot(n_chrom == length(chromosomes))

	# Translating the numeric values of the input file into chromosome names
	gapit_df$Chromosome <- chromosomes[gapit_df$Chromosome]

	# Coercing the gapit_df data.frame to a GRanges object
	output <- GenomicRanges::makeGRangesFromDataFrame(gapit_df,
							  start.field = "Position",
							  seqnames.field = "Chromosome",
							  end.field = "Position",
							  keep.extra.columns = TRUE)

	# Setting the seqlevels and seqlengths of the output GRanges from the .fai index
	GenomeInfoDb::seqlevels(output) <- GenomeInfoDb::seqlevels(fai_info)
	GenomeInfoDb::seqlengths(output) <- GenomeInfoDb::seqlengths(fai_info)

	# Pruning the levels according to the pattern
	if(!is.null(pattern)) GenomeInfoDb::seqlevels(output) <- grep(pattern, GenomeInfoDb::seqlevels(output), value = TRUE)

	# Adding a column for the -log10(corrected p-value) for downstream applications
	output$log10p <- -log10(gapit_df$P.value)

	# Returning the output
	return(output)
}

