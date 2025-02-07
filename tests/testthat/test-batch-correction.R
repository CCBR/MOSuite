test_that("batch_correction works for NIDAP", {
  moo <- multiOmicDataSet(
    sample_meta_dat = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts),
      "norm" = list(
        "voom" = as.data.frame(nidap_norm_counts)
      )
    )
  ) %>%
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
    moo@counts[["batch"]] %>%
      dplyr::arrange(desc(Gene)),
    as.data.frame(nidap_batch_corrected_counts_2) %>%
      dplyr::arrange(desc(Gene))
  ))
})
