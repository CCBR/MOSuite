test_that("plot_venn_diagram works with defaults", {
  expect_snapshot(
    p <- plot_venn_diagram(
      nidap_volcano_summary_dat,
      print_plots = FALSE,
      save_plots = TRUE
    )
  )
  # expect_equal(plot_venn_diagram(nidap_volcano_summary_dat), as.data.frame(nidap_venn_diagram_dat))
})
