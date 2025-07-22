test_that("plot_volcano_summary works on nidap dataset", {
  expect_snapshot(
    df_volc_sum <- plot_volcano_summary(nidap_deg_analysis, save_plots = TRUE)
  )
})
