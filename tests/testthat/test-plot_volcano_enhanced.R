test_that("plot_volcano_enhanced works on nidap dataset", {
  skip_on_ci()
  expect_snapshot(
    df_volc_enh <- plot_volcano_enhanced(
      nidap_deg_analysis,
      save_plots = FALSE,
      print_plots = FALSE
    )
  )
})

test_that("plot_volcano_enhanced returns a data frame", {
  expect_no_error(
    result <- plot_volcano_enhanced(
      nidap_deg_analysis,
      save_plots = FALSE,
      print_plots = FALSE
    )
  )

  expect_s3_class(result, "data.frame")
  expect_true(ncol(result) > 0)
  expect_true(nrow(result) > 0)
})

test_that("plot_volcano_enhanced respects num_features_to_label", {
  expect_no_error(
    result <- plot_volcano_enhanced(
      nidap_deg_analysis,
      num_features_to_label = 10,
      save_plots = FALSE,
      print_plots = FALSE
    )
  )

  expect_s3_class(result, "data.frame")
})

test_that("plot_volcano_enhanced forwards shared styling parameters", {
  options(mosuite_test_volcano_args = list())
  trace(
    EnhancedVolcano::EnhancedVolcano,
    tracer = quote(options(
      mosuite_test_volcano_args = append(
        getOption("mosuite_test_volcano_args"),
        list(list(
          labSize = labSize,
          labCol = labCol,
          col = col,
          cutoffLineCol = cutoffLineCol
        ))
      )
    )),
    print = FALSE
  )
  on.exit(untrace(EnhancedVolcano::EnhancedVolcano), add = TRUE)
  on.exit(options(mosuite_test_volcano_args = NULL), add = TRUE)

  custom_label <- tail(nidap_deg_analysis$Gene, 1)
  point_colors <- c("grey40", "orange", "dodgerblue", "firebrick")

  result <- plot_volcano_enhanced(
    nidap_deg_analysis,
    label_features = TRUE,
    custom_gene_list = custom_label,
    label_font_size = 7,
    default_label_color = "purple",
    custom_label_color = "darkgreen",
    color_of_signif_threshold_line = "cyan",
    color_of_non_significant_features = point_colors[1],
    color_of_logfold_change_threshold_line = point_colors[2],
    color_of_features_meeting_only_signif_threshold = point_colors[3],
    color_for_features_meeting_pvalue_and_foldchange_thresholds = point_colors[4],
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  captured_args <- getOption("mosuite_test_volcano_args")[[1]]
  expect_equal(captured_args$labSize, 7)
  expect_true(all(captured_args$labCol == "darkgreen"))
  expect_equal(captured_args$col, point_colors)
  expect_equal(captured_args$cutoffLineCol, "cyan")
})

test_that("plot_volcano_enhanced uses EnhancedVolcano default colors", {
  options(mosuite_test_volcano_args = list())
  trace(
    EnhancedVolcano::EnhancedVolcano,
    tracer = quote(options(
      mosuite_test_volcano_args = append(
        getOption("mosuite_test_volcano_args"),
        list(list(
          col = col,
          cutoffLineCol = cutoffLineCol
        ))
      )
    )),
    print = FALSE
  )
  on.exit(untrace(EnhancedVolcano::EnhancedVolcano), add = TRUE)
  on.exit(options(mosuite_test_volcano_args = NULL), add = TRUE)

  result <- plot_volcano_enhanced(
    nidap_deg_analysis,
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  captured_args <- getOption("mosuite_test_volcano_args")[[1]]
  expect_equal(captured_args$col, c("grey30", "forestgreen", "royalblue", "red2"))
  expect_equal(captured_args$cutoffLineCol, "black")
})

test_that("plot_volcano_enhanced displays selected genes", {
  options(mosuite_test_select_labels = list())
  trace(
    EnhancedVolcano::EnhancedVolcano,
    tracer = quote(options(
      mosuite_test_select_labels = append(
        getOption("mosuite_test_select_labels"),
        list(list(
          selectLab = selectLab,
          labCol = labCol
        ))
      )
    )),
    print = FALSE
  )
  on.exit(untrace(EnhancedVolcano::EnhancedVolcano), add = TRUE)
  on.exit(options(mosuite_test_select_labels = NULL), add = TRUE)

  selected_genes <- nidap_deg_analysis$Gene[1:2]
  result <- plot_volcano_enhanced(
    nidap_deg_analysis,
    label_features = TRUE,
    custom_gene_list = paste(selected_genes, collapse = ","),
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  captured_select_labels <- getOption("mosuite_test_select_labels")
  expect_true(length(captured_select_labels) > 0)
  expect_true(all(vapply(
    captured_select_labels,
    function(x) setequal(x$selectLab, selected_genes),
    logical(1)
  )))
  expect_true(all(vapply(
    captured_select_labels,
    function(x) all(x$labCol == "black"),
    logical(1)
  )))
})

test_that("plot_volcano_enhanced defaults invalid manual grid rows to one", {
  options(mosuite_test_grid_rows = list())
  trace(
    patchwork::wrap_plots,
    tracer = quote(options(
      mosuite_test_grid_rows = append(
        getOption("mosuite_test_grid_rows"),
        list(nrow)
      )
    )),
    print = FALSE
  )
  on.exit(untrace(patchwork::wrap_plots), add = TRUE)
  on.exit(options(mosuite_test_grid_rows = NULL), add = TRUE)

  invalid_rows <- list(0, -1, NA_real_, c(1, 2))
  for (invalid_row in invalid_rows) {
    options(mosuite_test_grid_rows = list())
    expect_no_error(
      result <- plot_volcano_enhanced(
        nidap_deg_analysis,
        use_default_grid_layout = FALSE,
        number_of_rows_in_grid_layout = invalid_row,
        save_plots = FALSE,
        print_plots = FALSE
      )
    )

    expect_s3_class(result, "data.frame")
    expect_equal(getOption("mosuite_test_grid_rows")[[1]], 1L)
  }
})

test_that("plot_volcano_enhanced scales grid output dimensions", {
  options(
    mosuite_test_grid_rows = list(),
    mosuite_test_plot_output_args = list()
  )
  trace(
    patchwork::wrap_plots,
    tracer = quote(options(
      mosuite_test_grid_rows = append(
        getOption("mosuite_test_grid_rows"),
        list(nrow)
      )
    )),
    print = FALSE
  )
  trace(
    ggplot2::ggsave,
    tracer = quote({
      options(
        mosuite_test_plot_output_args = append(
          getOption("mosuite_test_plot_output_args"),
          list(list(
            width = width,
            height = height,
            units = units,
            dpi = dpi,
            filename = filename
          ))
        )
      )
    }),
    print = FALSE
  )
  on.exit(untrace(patchwork::wrap_plots), add = TRUE)
  on.exit(untrace(ggplot2::ggsave), add = TRUE)
  on.exit(options(
    mosuite_test_grid_rows = NULL,
    mosuite_test_plot_output_args = NULL
  ), add = TRUE)

  plots_dir <- tempfile("volcano-grid-output-")
  dir.create(plots_dir)
  on.exit(unlink(plots_dir, recursive = TRUE), add = TRUE)

  cases <- list(
    list(nrows = 2, scale = TRUE, width = 100, height = 400),
    list(nrows = 1, scale = TRUE, width = 200, height = 200),
    list(nrows = 2, scale = FALSE, width = 100, height = 200)
  )

  for (case in cases) {
    options(
      mosuite_test_grid_rows = list(),
      mosuite_test_plot_output_args = list()
    )
    expect_no_error(
      result <- plot_volcano_enhanced(
        nidap_deg_analysis,
        use_default_grid_layout = FALSE,
        number_of_rows_in_grid_layout = case$nrows,
        scale_image_to_grid = case$scale,
        image_width = 100,
        image_height = 200,
        save_plots = TRUE,
        print_plots = FALSE,
        plots_subdir = plots_dir
      )
    )

    expect_s3_class(result, "data.frame")
    captured_grid_rows <- getOption("mosuite_test_grid_rows")[[1]]
    captured_output_args <- getOption("mosuite_test_plot_output_args")[[1]]
    expect_equal(captured_grid_rows, as.integer(case$nrows))
    expect_equal(captured_output_args$width, case$width)
    expect_equal(captured_output_args$height, case$height)
    expect_equal(captured_output_args$units, "px")
    expect_equal(captured_output_args$dpi, 300)
    expect_match(captured_output_args$filename, basename(plots_dir), fixed = TRUE)
    expect_match(captured_output_args$filename, "volcano_enhanced.png", fixed = TRUE)
  }
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

  expect_no_error(
    result <- plot_volcano_enhanced(
      moo,
      save_plots = FALSE,
      print_plots = FALSE
    )
  )

  expect_s3_class(result, "data.frame")
  expect_true(ncol(result) > 0)
  expect_true(nrow(result) > 0)
})
