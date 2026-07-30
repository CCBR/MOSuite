test_that("plot_volcano_summary works on nidap dataset", {
  expect_snapshot(
    df_volc_sum <- plot_volcano_summary(
      nidap_deg_analysis,
      save_plots = FALSE,
      print_plots = FALSE
    )
  )
  expect_s3_class(df_volc_sum, "data.frame")
  expect_true(nrow(df_volc_sum) > 0)
  expect_true(all(c("Gene", "Contrast") %in% colnames(df_volc_sum)))
})

test_that("plot_volcano_summary respects non-Gene feature ID column", {
  deg_analysis <- nidap_deg_analysis |>
    dplyr::rename(feature_id = Gene)

  df_volc_sum <- plot_volcano_summary(
    deg_analysis,
    feature_id_colname = "feature_id",
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(df_volc_sum, "data.frame")
  expect_true(nrow(df_volc_sum) > 0)
  expect_true("feature_id" %in% colnames(df_volc_sum))
  expect_false("Gene" %in% colnames(df_volc_sum))
})

test_that("plot_volcano_summary only forwards custom labels when requested", {
  options(mosuite_test_select_labels = list())
  trace(
    EnhancedVolcano::EnhancedVolcano,
    tracer = quote(options(
      mosuite_test_select_labels = append(
        getOption("mosuite_test_select_labels"),
        list(selectLab)
      )
    )),
    print = FALSE
  )
  on.exit(untrace(EnhancedVolcano::EnhancedVolcano), add = TRUE)
  on.exit(options(mosuite_test_select_labels = NULL), add = TRUE)

  custom_label <- tail(nidap_deg_analysis$Gene, 1)

  result <- plot_volcano_summary(
    nidap_deg_analysis,
    custom_gene_list = custom_label,
    add_features = FALSE,
    label_features = FALSE,
    num_features_to_label = 1,
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  captured_select_labels <- getOption("mosuite_test_select_labels")
  expect_false(custom_label %in% captured_select_labels[[1]])

  options(mosuite_test_select_labels = list())

  result <- plot_volcano_summary(
    nidap_deg_analysis,
    custom_gene_list = custom_label,
    add_features = TRUE,
    label_features = FALSE,
    num_features_to_label = 1,
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  captured_select_labels <- getOption("mosuite_test_select_labels")
  expect_true(custom_label %in% captured_select_labels[[1]])
})

test_that("plot_volcano_summary forwards shared styling parameters", {
  options(mosuite_test_volcano_args = list())
  trace(
    EnhancedVolcano::EnhancedVolcano,
    tracer = quote(options(
      mosuite_test_volcano_args = append(
        getOption("mosuite_test_volcano_args"),
        list(list(
          axisLabSize = axisLabSize,
          titleLabSize = titleLabSize,
          labSize = labSize,
          drawConnectors = drawConnectors,
          labCol = labCol,
          xlab = xlab,
          ylab = ylab,
          xlim = xlim,
          ylim = ylim,
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

  result <- plot_volcano_summary(
    nidap_deg_analysis,
    add_features = TRUE,
    custom_gene_list = custom_label,
    label_font_size = 7,
    draw_connectors = FALSE,
    axis_lab_size = 25,
    title_font_size = 26,
    change_lfc_name = "Custom fold change",
    change_sig_name = "Custom significance",
    use_custom_lab = TRUE,
    default_label_color = "purple",
    custom_label_color = "darkgreen",
    color_of_signif_threshold_line = "cyan",
    color_of_non_significant_features = point_colors[1],
    color_of_logfold_change_threshold_line = point_colors[2],
    color_of_features_meeting_only_signif_threshold = point_colors[3],
    color_for_features_meeting_pvalue_and_foldchange_thresholds = point_colors[
      4
    ],
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  captured_args <- getOption("mosuite_test_volcano_args")[[1]]
  expect_equal(captured_args$axisLabSize, 25)
  expect_equal(captured_args$titleLabSize, 26)
  expect_equal(captured_args$labSize, 7)
  expect_false(captured_args$drawConnectors)
  expect_equal(captured_args$xlab, "Custom fold change")
  expect_equal(captured_args$ylab, "Custom significance")
  expect_equal(captured_args$xlim, c(-8, 11))
  expect_equal(captured_args$ylim, c(0, 5))
  expect_true("darkgreen" %in% captured_args$labCol)
  expect_equal(captured_args$col, point_colors)
  expect_equal(captured_args$cutoffLineCol, "cyan")
})

test_that("plot_volcano_summary uses EnhancedVolcano default colors", {
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

  result <- plot_volcano_summary(
    nidap_deg_analysis,
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  captured_args <- getOption("mosuite_test_volcano_args")[[1]]
  expect_equal(
    captured_args$col,
    c("grey30", "forestgreen", "royalblue", "red2")
  )
  expect_equal(captured_args$cutoffLineCol, "black")
})

test_that("plot_volcano_summary uses comparison titles without subtitles", {
  options(mosuite_test_volcano_args = list())
  trace(
    EnhancedVolcano::EnhancedVolcano,
    tracer = quote(options(
      mosuite_test_volcano_args = append(
        getOption("mosuite_test_volcano_args"),
        list(list(
          title = title,
          subtitle = subtitle
        ))
      )
    )),
    print = FALSE
  )
  on.exit(untrace(EnhancedVolcano::EnhancedVolcano), add = TRUE)
  on.exit(options(mosuite_test_volcano_args = NULL), add = TRUE)

  result <- plot_volcano_summary(
    nidap_deg_analysis,
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  captured_args <- getOption("mosuite_test_volcano_args")
  expect_equal(
    vapply(captured_args, `[[`, character(1), "title"),
    c("B-A", "C-A", "B-C")
  )
  expect_true(all(vapply(
    captured_args,
    function(x) is.null(x$subtitle),
    logical(1)
  )))
})

test_that("plot_volcano_summary combines enhanced plots into one grid figure", {
  options(mosuite_test_plot_output_args = list())
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
  on.exit(untrace(ggplot2::ggsave), add = TRUE)
  on.exit(options(mosuite_test_plot_output_args = NULL), add = TRUE)

  plots_dir <- tempfile("volcano-summary-grid-output-")
  dir.create(plots_dir)
  on.exit(unlink(plots_dir, recursive = TRUE), add = TRUE)

  result <- plot_volcano_summary(
    nidap_deg_analysis,
    image_width = 1,
    image_height = 2,
    dpi = 100,
    draw_connectors = FALSE,
    use_default_grid_layout = FALSE,
    number_of_rows_in_grid_layout = 1,
    save_plots = TRUE,
    print_plots = FALSE,
    plots_subdir = plots_dir
  )

  expect_s3_class(result, "data.frame")
  captured_output_args <- getOption("mosuite_test_plot_output_args")
  expect_length(captured_output_args, 1)
  expect_equal(captured_output_args[[1]]$width, 300)
  expect_equal(captured_output_args[[1]]$height, 200)
  expect_equal(captured_output_args[[1]]$units, "px")
  expect_equal(captured_output_args[[1]]$dpi, 100)
  expect_match(
    captured_output_args[[1]]$filename,
    basename(plots_dir),
    fixed = TRUE
  )
  expect_match(
    captured_output_args[[1]]$filename,
    "volcano_summary.png",
    fixed = TRUE
  )
})

test_that("plot_volcano_summary displays selected genes", {
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
  result <- plot_volcano_summary(
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

test_that("plot_volcano_summary offsets labels when connectors are enabled", {
  options(mosuite_test_saved_summary_plot = list())
  trace(
    ggplot2::ggsave,
    tracer = quote(options(
      mosuite_test_saved_summary_plot = append(
        getOption("mosuite_test_saved_summary_plot"),
        list(list(plot = plot))
      )
    )),
    print = FALSE
  )
  on.exit(untrace(ggplot2::ggsave), add = TRUE)
  on.exit(options(mosuite_test_saved_summary_plot = NULL), add = TRUE)

  plots_dir <- tempfile("volcano-summary-label-offset-")
  dir.create(plots_dir)
  on.exit(unlink(plots_dir, recursive = TRUE), add = TRUE)

  volcano_data <- data.frame(
    Gene = c("gene_a", "gene_b", "not_selected"),
    `B-A_logFC` = c(2, -2, 0.1),
    `B-A_pval` = c(0.001, 0.002, 0.5),
    `B-A_tstat` = c(8, -7, 0.5),
    check.names = FALSE
  )

  result <- suppressWarnings(plot_volcano_summary(
    volcano_data,
    feature_id_colname = "Gene",
    signif_colname = "B-A_pval",
    label_features = TRUE,
    custom_gene_list = "gene_a,gene_b",
    draw_connectors = TRUE,
    add_deg_columns = "none",
    save_plots = TRUE,
    print_plots = FALSE,
    plots_subdir = plots_dir
  ))

  expect_s3_class(result, "data.frame")
  captured_plot <- getOption("mosuite_test_saved_summary_plot")[[1]]$plot
  layer_geoms <- vapply(
    captured_plot$layers,
    function(layer) class(layer$geom)[1],
    character(1)
  )
  expect_true("GeomTextRepel" %in% layer_geoms)
  expect_false("GeomText" %in% layer_geoms)
})

test_that("plot_volcano_summary works with multiOmicDataSet", {
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
  result <- plot_volcano_summary(
    moo,
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_true(ncol(result) > 0)
  expect_true(nrow(result) > 0)
  expect_true("Gene" %in% colnames(result))
  expect_true("Contrast" %in% colnames(result))
})

test_that("plot_volcano_summary auto-detects adjpval in preference to pval", {
  options(mosuite_test_volcano_args = list())
  trace(
    EnhancedVolcano::EnhancedVolcano,
    tracer = quote(options(
      mosuite_test_volcano_args = append(
        getOption("mosuite_test_volcano_args"),
        list(list(y = y))
      )
    )),
    print = FALSE
  )
  on.exit(untrace(EnhancedVolcano::EnhancedVolcano), add = TRUE)
  on.exit(options(mosuite_test_volcano_args = NULL), add = TRUE)

  # Dataset with both pval and adjpval columns
  df_both <- data.frame(
    Gene = letters[1:5],
    `B-A_logFC` = c(2, -2, 0.1, 1.5, -1.5),
    `B-A_pval` = c(0.001, 0.002, 0.5, 0.01, 0.02),
    `B-A_adjpval` = c(0.01, 0.02, 0.5, 0.05, 0.06),
    `B-A_tstat` = c(8, -7, 0.5, 5, -5),
    check.names = FALSE
  )

  result <- plot_volcano_summary(
    df_both,
    feature_id_colname = "Gene",
    add_deg_columns = "none",
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  captured_args <- getOption("mosuite_test_volcano_args")[[1]]
  expect_equal(captured_args$y, "B-A_adjpval")
})

test_that("plot_volcano_summary change_colname restricts contrasts plotted", {
  result <- plot_volcano_summary(
    nidap_deg_analysis,
    change_colname = "B-A_logFC",
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(result$Contrast == "B-A"))
})

test_that("plot_volcano_summary explicit signif_colname as full column name is respected", {
  options(mosuite_test_volcano_args = list())
  trace(
    EnhancedVolcano::EnhancedVolcano,
    tracer = quote(options(
      mosuite_test_volcano_args = append(
        getOption("mosuite_test_volcano_args"),
        list(list(y = y))
      )
    )),
    print = FALSE
  )
  on.exit(untrace(EnhancedVolcano::EnhancedVolcano), add = TRUE)
  on.exit(options(mosuite_test_volcano_args = NULL), add = TRUE)

  result <- plot_volcano_summary(
    nidap_deg_analysis,
    change_colname = "B-A_logFC",
    signif_colname = "B-A_pval",
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  captured_args <- getOption("mosuite_test_volcano_args")[[1]]
  expect_equal(captured_args$y, "B-A_pval")
})

test_that("plot_volcano_summary errors when no significance columns found", {
  bad_df <- data.frame(
    Gene = letters[1:3],
    `B-A_logFC` = c(1, -1, 0.1),
    `B-A_tstat` = c(5, -4, 0.2),
    check.names = FALSE
  )
  expect_error(
    plot_volcano_summary(bad_df, save_plots = FALSE, print_plots = FALSE),
    regexp = "auto-detect"
  )
})
