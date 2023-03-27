# Overview

The `gwask` package includes functions for performing post-GWAS analysis
of [*k*-mer GWAS](https://github.com/voichek/kmersGWAS).

# Installation

## Dependencies

Some Bioconductor packages will need to be installed for this package to work:

* Biostrings
* GenomeInfoDb
* GenomicRanges
* IRanges
* MatrixGenerics
* Rsamtools
* S4Vectors
* VariantAnnotation

[GAPIT3](https://github.com/jiabowang/GAPIT) also needs to be installed.

Other required packages should be pulled in automatically from CRAN.

## Installing the package

Package sources can be downloaded from GitHub by running `git clone https://github.com/malemay/gwask`.
Then, running the following command in `R` should install the package from source:

	install.packages("gwask", repos = NULL, type = "source")

# Documentation

There is no vignette availble at the moment for `gwask`, however all functions
included in the package are individually documented. A complete list of the
available functions can be found here:

* `adjust_gaps`: Adjust plotting position according to alignment gaps
* `cluster_haplotypes`: Cluster the most similar haplotypes together
* `cluster_ld`: Greedy clustering of an LD matrix
* `extract_signals`: Extract signal ranges from GWAS results
* `fill_gaps`: Fill the deleted positions with dashes
* `format_gapit_gwas`: Format GAPIT GWAS results for generating Manhattan plots
* `format_haplotypes`: Format haplotypes for plotting using k-mer overlap information
* `format_kmer_gwas`: Format k-mer GWAS results for downstream analyses
* `gapit_vcf`: Launch a MLM GWAS analysis on a VCF file read with VariantAnnotation::readVcf using GAPIT
* `get_haplotypes`: Extract and filter the haplotypes from a set of sequences
* `gg_color_hue`: Select colors for a discrete scale as in ggplot2
* `gg_hue`: Generating a vector of colors comparable to those used by ggplot2
* `grid.colorscale`: Plot the color scale used in a haplotype plot
* `grid.haplotypes`: Plot a set of haplotypes using grid functions
* `grid.phenotable`: Plot a contingency table of observed phenotypes and haplotypes
* `is_valid_ld`: Check the validity of an LD matrix
* `kmer_ld`: Compute the pairwise LD between k-mers based on their presence/absence
* `ld_plot`: Plot an LD matrix using grid functions
* `ld_sort`: Re-arrange the samples in an LD matrix according to user-specified criteria
* `link_phenotypes`: Link haplotypes to their observed phenotypes
* `mafft_align`: Align a set of sequences using mafft
* `manhattanGrob`: A function that returns a graphical object (grob) representing a manhattan plot
* `map_color`: Map a set of numeric values onto a color palette
* `match_kmers`: Match a set of k-mers to positions on sequences
* `nucdiff`: Find positions that differ between two sequences
* `pvalueGrob`: A grob representing p-values of a GWAS analysis at a locus to plot using grid functions
* `pvalue_tx_grob`: Arrange a transcript grob and several p-value grobs in the same plot
* `read_fasta`: Read a fasta file into a named character vector
* `read_kmer_pvalues`: Read the (sorted) k-mer p-values from the output of a k-mer GWAS analysis
* `subsample_kmers`: Subsample from a set of significant k-mers prior to computing LD
* `transcriptGrob`: Generate a grob representing a transcript using grid functions
* `transcriptsGrob`: Generate a grob of all the possible transcripts for genes in a genomic region using grid
* `vcf_to_gapit`: Convert VCF records read with VariantAnnotation::readVCF into a format usable by GAPIT

# Citation

The reference to the `gwask` package will be uploaded shortly.

