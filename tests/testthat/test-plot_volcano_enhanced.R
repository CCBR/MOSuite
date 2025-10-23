test_that("plot_volcano_enhanced works on nidap dataset", {
  expect_snapshot(
    df_volc_enh <- plot_volcano_enhanced(
      nidap_deg_analysis,
      save_plots = FALSE,
      print_plots = FALSE
    )
  )
})
