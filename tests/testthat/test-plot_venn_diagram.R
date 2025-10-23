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
test_that("plot_venn_diagram raises condition for empty df", {
  expect_error(plot_venn_diagram(structure(
    list(
      GeneName = character(0),
      Contrast = character(0),
      FC = numeric(0),
      logFC = numeric(0),
      tstat = numeric(0),
      pval = numeric(0),
      adjpval = numeric(0)
    ),
    class = "data.frame",
    row.names = integer(0)
  )), "Dataframe is empty")
})
