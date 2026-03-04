test_that("plot_volcano_enhanced works on nidap dataset", {
  expect_snapshot(
    df_volc_enh <- plot_volcano_enhanced(
      nidap_deg_analysis,
      save_plots = FALSE,
      print_plots = FALSE
    )
  )
})

test_that("plot_volcano_enhanced returns a data frame", {
  result <- plot_volcano_enhanced(
    nidap_deg_analysis,
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_true(ncol(result) > 0)
  expect_true(nrow(result) > 0)
})

test_that("plot_volcano_enhanced respects num_features_to_label", {
  result <- plot_volcano_enhanced(
    nidap_deg_analysis,
    num_features_to_label = 10,
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
})

test_that("plot_volcano_enhanced works with multiOmicDataSet", {
  # Create a multiOmicDataSet with differential analysis results
  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts)
    ),
    analyses_lst = list(
      diff = nidap_deg_analysis_2
    )
  )

  # Test that it returns a data frame
  result <- plot_volcano_enhanced(
    moo,
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_true(ncol(result) > 0)
  expect_true(nrow(result) > 0)
})
