test_that("differential analysis works for NIDAP", {
  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts)
    )
  ) %>%
    analyze_diff_counts(
      count_type = "filt",
      sub_count_type = NULL,
      covariates_colnames = c("Group", "Batch"),
      contrast_colname = c("Group"),
      contrasts = c("B-A", "C-A", "B-C"),
      voom_normalization_method = "quantile",
    )

  x <- nidap_deg_analysis %>%
    dplyr::arrange(Gene) %>%
    dplyr::select(order(colnames(.))) %>%
    as.data.frame()

  y <- moo@analyses$diff %>%
    dplyr::arrange(Gene) %>%
    dplyr::select(order(colnames(.)))

  equal_dfs(x, y)
})
