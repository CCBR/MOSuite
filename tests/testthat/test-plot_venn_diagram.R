test_that("plot_venn_diagram works with defaults", {
  expect_snapshot(
    p <- plot_venn_diagram(
      nidap_volcano_summary_dat,
      print_plots = FALSE,
      save_plots = TRUE
    )
  )
  expect_equal(
    plot_venn_diagram(nidap_volcano_summary_dat),
    as.data.frame(nidap_venn_diagram_dat)
  )
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

test_that("intersection matrix assignment avoids recursive evaluation error", {
  # This test demonstrates the fix for the recursive default argument reference error
  # The error occurred with this pattern: Intersection <- sapply(colnames(Intersection), function(x) Intersection[, x])
  # The fix uses a temporary variable: intersection_matrix <- Intersection; Intersection <- sapply(colnames(intersection_matrix), ...)

  # Call plot_venn_diagram directly to ensure the fix works in practice
  expect_no_error({
    result <- plot_venn_diagram(
      nidap_volcano_summary_dat,
      print_plots = FALSE,
      save_plots = FALSE
    )
  })

  # Verify the result has the expected structure
  expect_s3_class(result, "data.frame")
  expect_true("Gene" %in% colnames(result))
  expect_true("Intersection" %in% colnames(result))
  expect_true("Id" %in% colnames(result))
  expect_true("Size" %in% colnames(result))
})
