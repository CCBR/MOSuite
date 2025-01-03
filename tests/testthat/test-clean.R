test_that("clean_raw_counts works", {
  moo_nidap <- create_multiOmicDataSet_from_dataframes(
    sample_meta_dat = as.data.frame(nidap_sample_metadata),
    count_dat = as.data.frame(nidap_raw_counts),
    sample_id_colname = "Sample"
  ) %>%
    clean_raw_counts(
      sample_names_column = "Sample",
      gene_names_column = "GeneName"
    )
  expect_true(equal_dfs(moo_nidap@counts[["clean"]], as.data.frame(nidap_clean_raw_counts)))
})

test_that("check_sample_names works for dataframes & tibbles", {
  expect_true(check_sample_names(nidap_raw_counts, nidap_sample_metadata, sample_names_column = "Sample"))
  expect_true(check_sample_names(as.data.frame(nidap_raw_counts), as.data.frame(nidap_sample_metadata), sample_names_column = "Sample"))
  expect_error(
    check_sample_names(nidap_raw_counts %>% dplyr::select(-A1),
      nidap_sample_metadata,
      sample_names_column = "Sample"
    ),
    regexp = "The following sample names are in the metadata but not the raw counts"
  )
})
