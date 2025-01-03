test_that("normalize works", {
  moo <- multiOmicDataSet(
    sample_meta_dat = as.data.frame(nidap_sample_metadata),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts)
    )
  ) %>%
    normalize(
      gene_names_column = "Gene",
      columns_to_include = c("Gene", "A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"),
      sample_names_column = "Sample",
      group_column = "Group",
      label_column = "Label"
    )
  expect_true(equal_dfs(
    moo@counts[["norm"]][["voom"]] %>%
      dplyr::arrange(desc(Gene)),
    as.data.frame(nidap_norm_counts) %>%
      dplyr::arrange(desc(Gene))
  ))
})
