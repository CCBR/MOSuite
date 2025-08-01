#' RSEM expected gene counts
#'
#' @format ## `gene_counts`
#' A data frame with columns 'gene_id', 'GeneName', and a column for each sample's expected count.
#' @keywords data
#' @source Generated by running RENEE v2.5.8 on the
#' [test dataset](https://github.com/CCBR/RENEE/tree/e08f7db6c6e638cfd330caa182f64665d2ef37fa/.tests)
"gene_counts"

#' Sample metadata for the NIDAP test dataset
#' @keywords data
"nidap_sample_metadata"

#' Raw counts for the NIDAP test dataset
#' Pairs with `nidap_sample_metadata`.
#' @keywords data
"nidap_raw_counts"

#' Clean raw counts for the NIDAP test dataset.
#' The result of running `clean_raw_counts()` on `nidap_raw_counts`.
#' @keywords data
"nidap_clean_raw_counts"

#' Filtered counts for the NIDAP test dataset.
#' The result of running `filter_counts()` on `nidap_clean_raw_counts`.
#' @keywords data
"nidap_filtered_counts"

#' Normalized counts for the NIDAP test dataset.
#' The result of running `normalize_counts()` on `nidap_filtered_counts`.
#' @keywords data
"nidap_norm_counts"

#' Batch-corrected counts for the NIDAP test dataset.
#' @keywords data
"nidap_batch_corrected_counts"

#' Batch-corrected counts for the NIDAP test dataset.
#' The result of running `batch_correct_counts()` on `nidap_norm_counts`.
#' @keywords data
"nidap_batch_corrected_counts_2"


#' Differential gene expression analysis for the NIDAP test dataset.
#' @keywords data
"nidap_deg_analysis"


#' Differential gene expression analysis for the NIDAP test dataset.
#' The result of running `diff_counts()` on `nidap_filtered_counts`.
#' @keywords data
"nidap_deg_analysis_2"

#' List of differentially expressed genes from the NIDAP test dataset using
#' default parameters with `filter_diff()`.
#' @keywords data
"nidap_deg_gene_list"

#' Summarized differential expression analysis for input to venn diagram
#' @keywords data
"nidap_volcano_summary_dat"

#' Output data from venn diagram.
#' The result of running `plot_venn_diagram()` on `nidap_volcano_summary_dat`
#' @keywords data
"nidap_venn_diagram_dat"
