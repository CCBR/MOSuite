test_that("clean_raw_counts works for NIDAP data", {
  moo_nidap <- create_multiOmicDataSet_from_dataframes(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    counts_dat = as.data.frame(nidap_raw_counts)
  ) %>%
    clean_raw_counts(print_plots = TRUE)
  expect_true(equal_dfs(
    moo_nidap@counts[["clean"]] %>%
      dplyr::rename(Gene = GeneName),
    as.data.frame(nidap_clean_raw_counts)
  ))
})

test_that("clean_raw_counts works for RENEE data", {
  moo <- create_multiOmicDataSet_from_dataframes(
    readr::read_tsv(
      system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite")
    ),
    gene_counts
  ) %>%
    clean_raw_counts()
  expect_equal(
    head(moo@counts$clean),
    structure(
      list(
        gene_id = c(
          "ENSG00000121410.11",
          "ENSG00000268895.5",
          "ENSG00000148584.15",
          "ENSG00000175899.14",
          "ENSG00000245105.3",
          "ENSG00000166535.20"
        ),
        KO_S3 = c(0, 0, 0, 0, 0, 0),
        KO_S4 = c(0, 0, 0, 0, 0, 0),
        WT_S1 = c(0, 0, 0, 0, 0, 0),
        WT_S2 = c(0, 0, 0, 0, 0, 0)
      ),
      row.names = c(NA, 6L),
      class = "data.frame"
    )
  )
  expect_equal(
    tail(moo@counts$clean),
    structure(
      list(
        gene_id = c(
          "ENSG00000232242.2",
          "ENSG00000162378.13",
          "ENSG00000159840.16",
          "ENSG00000274572.1",
          "ENSG00000074755.15",
          "ENSG00000272920.1"
        ),
        KO_S3 = c(0, 0, 0, 0, 0, 0),
        KO_S4 = c(0, 0, 0, 0, 0, 0),
        WT_S1 = c(0, 0, 0, 0, 0, 0),
        WT_S2 = c(0, 0, 0, 0, 0, 0)
      ),
      row.names = 58924:58929,
      class = "data.frame"
    )
  )
})
