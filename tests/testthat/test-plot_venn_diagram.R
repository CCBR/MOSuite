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
  expect_error(
    plot_venn_diagram(structure(
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
    )),
    "Dataframe is empty"
  )
})

test_that("plot_venn_diagram handles intersection matrix without recursive evaluation error", {
  # This test ensures that the fix for the recursive default argument reference error works
  # The error occurred when Intersection was being reassigned while still being referenced
  expect_no_error({
    result <- plot_venn_diagram(
      nidap_volcano_summary_dat,
      print_plots = FALSE,
      save_plots = FALSE
    )
  })
  # Verify the result is a data frame with expected columns
  expect_s3_class(result, "data.frame")
  expect_true("Gene" %in% colnames(result))
  expect_true("Intersection" %in% colnames(result))
  expect_true("Id" %in% colnames(result))
  expect_true("Size" %in% colnames(result))
})
