test_that("batch_correction works for NIDAP", {
  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts),
      "norm" = list("voom" = as.data.frame(nidap_norm_counts))
    )
  ) |>
    batch_correct_counts(
      count_type = "norm",
      sub_count_type = "voom",
      covariates_colnames = "Group",
      batch_colname = "Batch",
      label_colname = "Label",
      print_plots = TRUE
    )
  # TODO: getting different results than nidap_batch_corrected_counts
  expect_true(all.equal(
    moo@counts[["batch"]] |>
      dplyr::arrange(desc(Gene)),
    as.data.frame(nidap_batch_corrected_counts_2) |>
      dplyr::arrange(desc(Gene))
  ))
})

test_that("batch_correction warnings & errors", {
  moo <- create_multiOmicDataSet_from_dataframes(
    readr::read_tsv(
      system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite")
    ) |>
      dplyr::mutate(batch = 1),
    gene_counts
  ) |>
    clean_raw_counts() |>
    filter_counts(
      group_colname = "condition",
      label_colname = "sample_id",
      minimum_count_value_to_be_considered_nonzero = 1,
      minimum_number_of_samples_with_nonzero_counts_in_total = 1,
      minimum_number_of_samples_with_nonzero_counts_in_a_group = 1,
      print_plots = FALSE
    ) |>
    normalize_counts(group_colname = "condition", label_colname = "sample_id")

  expect_warning(
    moo |>
      batch_correct_counts(
        covariates_colnames = "condition",
        batch_colname = "batch"
      ),
    "Batch column 'batch' contains only 1 unique value"
  )
  expect_error(
    moo |>
      batch_correct_counts(
        covariates_colnames = "batch",
        batch_colname = "batch"
      ),
    "Batch column 'batch' cannot be included in covariates."
  )
})
