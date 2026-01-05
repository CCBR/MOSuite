testthat::test_that("aggregate_duplicate_gene_names returns collapsed dfout", {
  # Minimal counts table with duplicated feature IDs in the first column
  counts_dat <- data.frame(
    gene_id = c("A", "A", "B"),
    sample1 = c(1, 2, 3),
    sample2 = c(4, 5, 6),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  out <- MOSuite:::aggregate_duplicate_gene_names(
    counts_dat = counts_dat,
    gene_name_column_to_use_for_collapsing_duplicates = "gene_id",
    aggregate_rows_with_duplicate_gene_names = TRUE,
    split_gene_name = FALSE
  )
  
  # Should collapse A+A into one row, leaving A and B
  testthat::expect_equal(nrow(out), 2)
  testthat::expect_equal(sum(duplicated(out$gene_id)), 0)
  
  # Counts should be summed across duplicates
  a_row <- out[out$gene_id == "A", , drop = FALSE]
  testthat::expect_equal(a_row$sample1, 3)
  testthat::expect_equal(a_row$sample2, 9)
})
