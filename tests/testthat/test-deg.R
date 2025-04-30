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
    diff_counts(
      count_type = "filt",
      sub_count_type = NULL,
      sample_id_colname = "Sample",
      feature_id_colname = "Gene",
      columns_to_include = c("Gene", "A1", "A2", "A3", "B1", "B2", "B3", "C1", "C3"),
      covariates_colnames = c("Group", "Batch"),
      contrast_colname = c("Group"),
      contrasts = c("B-A", "C-A", "B-C"),
      voom_normalization_method = "quantile",
    )

  # x <- nidap_deg_analysis %>%
  #   dplyr::arrange(Gene) %>%
  #   dplyr::select(order(colnames(.))) %>% dplyr::select(-C2) %>%
  #   as.data.frame()
  # y <- moo@analyses$diff %>%
  #   dplyr::arrange(Gene) %>%
  #   dplyr::select(order(colnames(.)))
  # equal_dfs(x, y)

  expect_equal(
    moo@analyses$diff %>%
      dplyr::arrange(Gene) %>%
      dplyr::select(order(colnames(.))),
    nidap_deg_analysis_2 %>%
      dplyr::arrange(Gene) %>%
      dplyr::select(order(colnames(.)))
  )
})
