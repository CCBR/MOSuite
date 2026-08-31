make_venn_test_moo <- function() {
  contrast_ba <- data.frame(
    Gene = c("both", "ba_only", "p_boundary", "fc_boundary"),
    logFC = c(2, 2, 2, 1),
    adjpval = c(0.01, 0.01, 0.05, 0.01)
  )
  contrast_ca <- data.frame(
    Gene = c("both", "ca_only", "p_boundary", "fc_boundary"),
    logFC = c(-2, -2, -2, -1),
    adjpval = c(0.01, 0.01, 0.05, 0.01)
  )
  return(
    multiOmicDataSet(
      sample_metadata = data.frame(Sample = "S1"),
      anno_dat = data.frame(),
      counts_lst = list(raw = data.frame(Gene = "seed", S1 = 1)),
      analyses_lst = list(diff = list("B-A" = contrast_ba, "C-A" = contrast_ca))
    )
  )
}

test_that("plot_venn_diagram builds intersections directly from DEG MOO results", {
  result <- plot_venn_diagram(
    make_venn_test_moo(),
    print_plots = FALSE,
    save_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_setequal(result$Gene, c("both", "ba_only", "ca_only"))
  expect_false(any(result$Gene %in% c("p_boundary", "fc_boundary")))
  expect_true("(B-A ∩ C-A)" %in% result$Intersection)
})

test_that("plot_venn_diagram respects selected MOO contrasts", {
  result <- plot_venn_diagram(
    make_venn_test_moo(),
    select_contrasts = "B-A",
    print_plots = FALSE,
    save_plots = FALSE
  )

  expect_setequal(result$Gene, c("both", "ba_only"))
  expect_true(all(result$Intersection == "(B-A)"))
})

test_that("plot_venn_diagram errors for missing contrast DEG columns", {
  expect_error(
    plot_venn_diagram(
      make_venn_test_moo(),
      signif_colname = "missing_column",
      print_plots = FALSE,
      save_plots = FALSE
    ),
    "Required DEG columns not found"
  )
})

test_that("plot_venn_diagram uses an Intersection plot for more than five contrasts", {
  moo <- make_venn_test_moo()
  moo@analyses$diff <- rep(moo@analyses$diff["B-A"], 6)
  names(moo@analyses$diff) <- paste0("contrast_", seq_len(6))
  for (contrast_name in names(moo@analyses$diff)) {
    moo@analyses$diff[[contrast_name]] <- data.frame(
      Gene = c("both", contrast_name),
      logFC = c(2, 2),
      adjpval = c(0.01, 0.01)
    )
  }

  expect_no_error(
    suppressWarnings(
      plot_venn_diagram(
        moo,
        plot_type = "Venn diagram",
        print_plots = FALSE,
        save_plots = FALSE
      )
    )
  )
})
