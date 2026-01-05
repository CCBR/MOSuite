# training data set from the NIDAP Bulk RNA-seq workflow
moo <- multiOmicDataSet(
  sample_metadata = as.data.frame(nidap_sample_metadata),
  anno_dat = data.frame(),
  counts_lst = list(
    "raw" = as.data.frame(nidap_raw_counts),
    "clean" = as.data.frame(nidap_clean_raw_counts),
    "filt" = as.data.frame(nidap_filtered_counts),
    "norm" = list("voom" = as.data.frame(nidap_norm_counts))
  )
) %>%
  batch_correct_counts(
    count_type = "norm",
    sub_count_type = "voom",
    covariates_colname = "Group",
    batch_colname = "Batch",
    label_colname = "Label"
  ) %>%
  diff_counts(
    count_type = "filt",
    sub_count_type = NULL,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    covariates_colnames = c("Group", "Batch"),
    contrast_colname = c("Group"),
    contrasts = c("B-A", "C-A", "B-C"),
    voom_normalization_method = "quantile",
  ) %>%
  filter_diff(
    significance_column = "adjpval",
    significance_cutoff = 0.05,
    change_column = "logFC",
    change_cutoff = 1,
    filtering_mode = "any",
    include_estimates = c("FC", "logFC", "tstat", "pval", "adjpval"),
    round_estimates = TRUE,
    rounding_decimal_for_percent_cells = 0,
    contrast_filter = "none",
    contrasts = c(),
    groups = c(),
    groups_filter = "none",
    plot_type = "bar",
    save_plots = FALSE,
    print_plots = FALSE
  )


nidap_sample_metadata <- readr::read_csv(
  system.file(
    "extdata",
    "nidap",
    "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz",
    package = "MOSuite"
  )
)
usethis::use_data(nidap_sample_metadata, overwrite = TRUE)

nidap_raw_counts <- readr::read_csv(system.file(
  "extdata",
  "nidap",
  "Raw_Counts.csv.gz",
  package = "MOSuite"
))
usethis::use_data(nidap_raw_counts, overwrite = TRUE)

nidap_clean_raw_counts <- readr::read_csv(system.file(
  "extdata",
  "nidap",
  "Clean_Raw_Counts.csv.gz",
  package = "MOSuite"
))
usethis::use_data(nidap_clean_raw_counts, overwrite = TRUE)

nidap_filtered_counts <- readr::read_csv(system.file(
  "extdata",
  "nidap",
  "Filtered_Counts.csv.gz",
  package = "MOSuite"
))
usethis::use_data(nidap_filtered_counts, overwrite = TRUE)

nidap_norm_counts <- readr::read_csv(system.file(
  "extdata",
  "nidap",
  "Normalized_Counts.csv.gz",
  package = "MOSuite"
))
usethis::use_data(nidap_norm_counts, overwrite = TRUE)

nidap_batch_corrected_counts <- readr::read_csv(
  system.file(
    "extdata",
    "nidap",
    "Batch_Corrected_Counts.csv.gz",
    package = "MOSuite"
  )
)
usethis::use_data(nidap_batch_corrected_counts, overwrite = TRUE)

nidap_batch_corrected_counts_2 <- moo@counts[["batch"]]
usethis::use_data(nidap_batch_corrected_counts_2, overwrite = TRUE)

nidap_deg_analysis <- readr::read_csv(system.file(
  "extdata",
  "nidap",
  "DEG_Analysis.csv.gz",
  package = "MOSuite"
))
usethis::use_data(nidap_deg_analysis, overwrite = TRUE)


nidap_deg_analysis_2 <- moo@analyses$diff
usethis::use_data(nidap_deg_analysis_2, overwrite = TRUE)

# nidap_deg_gene_list <- readr::read_csv(system.file("extdata", "nidap", "DEG_Gene_List.csv.gz", package = "MOSuite"))
# usethis::use_data(nidap_deg_gene_list, overwrite = TRUE)

nidap_deg_gene_list <- moo@analyses$diff_filt
usethis::use_data(nidap_deg_gene_list, overwrite = TRUE)

nidap_volcano_summary_dat <- readr::read_csv(system.file(
  "extdata",
  "nidap",
  "Volcano_Summary.csv.gz",
  package = "MOSuite"
))
usethis::use_data(nidap_volcano_summary_dat, overwrite = TRUE)

nidap_venn_diagram_dat <- readr::read_csv(system.file(
  "extdata",
  "nidap",
  "Venn_Diagram.csv.gz",
  package = "MOSuite"
))
usethis::use_data(nidap_venn_diagram_dat, overwrite = TRUE)

## re-compress Rda files
# tools::resaveRdaFiles('data')
