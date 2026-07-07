test_that("plot_read_depth works on moo & dataframes", {
  moo <- multiOmicDataSet(
    sample_metadata = nidap_sample_metadata,
    anno_dat = data.frame(),
    counts_lst = list("raw" = nidap_raw_counts)
  )
  expect_equal(
    plot_read_depth(moo, "raw"),
    plot_read_depth(nidap_raw_counts)
  )
})

test_that("plot_read_depth accepts extra args via moo dispatch without error", {
  moo <- multiOmicDataSet(
    sample_metadata = nidap_sample_metadata,
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = nidap_raw_counts,
      "clean" = nidap_clean_raw_counts
    )
  )
  expect_no_error(
    plot_read_depth(
      moo,
      count_type = "clean",
      sample_id_colname = "Sample",
      feature_id_colname = "Gene"
    )
  )
})
