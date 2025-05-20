test_that("E2E workflow succeeds for RENEE data", {
  gene_counts_tsv <- system.file("extdata",
    "RSEM.genes.expected_count.all_samples.txt.gz",
    package = "MOSuite"
  )
  metadata_tsv <- system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite")

  expect_snapshot(
    moo <- create_multiOmicDataSet_from_files(
      sample_meta_filepath = metadata_tsv,
      feature_counts_filepath = gene_counts_tsv
    ) %>%
      clean_raw_counts() %>%
      filter_counts(
        group_colname = "condition",
        label_colname = "sample_id",
        minimum_count_value_to_be_considered_nonzero = 1,
        minimum_number_of_samples_with_nonzero_counts_in_total = 1,
        minimum_number_of_samples_with_nonzero_counts_in_a_group = 1,
      ) %>%
      normalize_counts(group_colname = "condition", label_colname = "sample_id") %>%
      diff_counts(
        covariates_colnames = "condition",
        contrast_colname = "condition",
        contrasts = c("knockout-wildtype")
      )
  )
})

test_that("E2E workflow succeeds for NIDAP data", {
  expect_snapshot(moo_nidap <- create_multiOmicDataSet_from_dataframes(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    counts_dat = as.data.frame(nidap_raw_counts)
  ) %>%
    clean_raw_counts() %>%
    filter_counts(group_colname = "Group") %>%
    normalize_counts(group_colname = "Group") %>%
    batch_correct_counts(
      covariates_colname = "Group",
      batch_colname = "Batch",
      label_colname = "Label"
    ) %>%
    diff_counts(
      count_type = "filt",
      sub_count_type = NULL,
      sample_id_colname = "Sample",
      feature_id_colname = "GeneName",
      covariates_colnames = c("Group", "Batch"),
      contrast_colname = c("Group"),
      contrasts = c("B-A", "C-A", "B-C"),
      input_in_log_counts = FALSE,
      return_mean_and_sd = FALSE,
      return_normalized_counts = TRUE,
      voom_normalization_method = "quantile",
    ))
})
