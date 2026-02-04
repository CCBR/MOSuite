test_that("clean_raw_counts works for NIDAP data", {
  moo_nidap <- create_multiOmicDataSet_from_dataframes(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    counts_dat = as.data.frame(nidap_raw_counts)
  ) %>%
    clean_raw_counts(print_plots = TRUE)

  actual <- moo_nidap@counts[["clean"]] %>%
    dplyr::rename(Gene = GeneName) %>%
    as.data.frame()

  expected <- as.data.frame(nidap_clean_raw_counts)

  cmp <- all.equal(actual, expected, check.attributes = FALSE)
  expect_true(isTRUE(cmp), info = paste(cmp, collapse = "\n"))
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

test_that("aggregate_duplicate_gene_names returns collapsed dfout", {
  counts_dat <- data.frame(
    gene_id = c("A", "A", "B"),
    sample1 = c(1, 2, 3),
    sample2 = c(4, 5, 6),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Case 1: aggregation enabled
  out <- MOSuite:::aggregate_duplicate_gene_names(
    counts_dat = counts_dat,
    gene_name_column_to_use_for_collapsing_duplicates = "gene_id",
    aggregate_rows_with_duplicate_gene_names = TRUE,
    split_gene_name = FALSE
  )

  expect_equal(nrow(out), 2)
  expect_equal(sum(duplicated(out$gene_id)), 0)

  a_row <- out[out$gene_id == "A", , drop = FALSE]
  expect_equal(a_row$sample1, 3)
  expect_equal(a_row$sample2, 9)

  # Case 2: aggregation disabled
  out_noagg <- MOSuite:::aggregate_duplicate_gene_names(
    counts_dat = counts_dat,
    gene_name_column_to_use_for_collapsing_duplicates = "gene_id",
    aggregate_rows_with_duplicate_gene_names = FALSE,
    split_gene_name = FALSE
  )

  expect_equal(nrow(out_noagg), 3)
  expect_equal(sum(duplicated(out_noagg$gene_id)), 1)
  expect_equal(out_noagg$sample1, counts_dat$sample1)
  expect_equal(out_noagg$sample2, counts_dat$sample2)
})
