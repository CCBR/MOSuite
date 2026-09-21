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
    color_for_features_meeting_pvalue_and_foldchange_thresholds = point_colors[
      4
    ],
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

test_that("build_volcano_plot_data colors exact threshold boundaries as non-significant", {
  diff_dat <- data.frame(
    Gene = c("both_boundary", "p_boundary", "fc_boundary", "neither"),
    annotation = c("a", "b", "c", "d"),
    `B-A_logFC` = c(1, 0.5, 1, 0.5),
    `B-A_pval` = c(0.05, 0.05, 0.06, 0.06),
    check.names = FALSE
  )

  volcano_data <- build_volcano_plot_data(
    diff_dat = diff_dat,
    label_col = "Gene",
    change_colname = "B-A_logFC",
    signif_colname = "B-A_pval",
    contrast_idx = 1,
    use_custom_lab = FALSE,
    change_lfc_name = "log2FC",
    change_sig_name = "p-value",
    value_to_sort_the_output_dataset = "p-value",
    label_significant_features_only = TRUE,
    signif_threshold = 0.05,
    change_threshold = 1,
    num_features_to_label = 30,
    custom_gene_list = "",
    label_features = FALSE,
    default_label_color = "black",
    custom_label_color = "black",
    color_of_non_significant_features = "grey30",
    color_of_logfold_change_threshold_line = "forestgreen",
    color_of_features_meeting_only_signif_threshold = "royalblue",
    color_for_features_meeting_pvalue_and_foldchange_thresholds = "red2"
  )

  point_colors_by_class <- volcano_data$custom_colors
  expect_equal(
    unname(point_colors_by_class[
      names(point_colors_by_class) == "Not significant"
    ]),
    rep("grey30", 4)
  )
  expect_true(
    all(
      point_colors_by_class[
        names(point_colors_by_class) == "Significant only"
      ] ==
        "royalblue"
    )
  )
  expect_true(
    all(
      point_colors_by_class[
        names(point_colors_by_class) == "Fold change only"
      ] ==
        "forestgreen"
    )
  )
})

test_that("plot_volcano_enhanced labels only strict threshold hits", {
  options(mosuite_test_volcano_label_args = list())
  trace(
    EnhancedVolcano::EnhancedVolcano,
    tracer = quote(options(
      mosuite_test_volcano_label_args = append(
        getOption("mosuite_test_volcano_label_args"),
        list(list(selectLab = selectLab))
      )
    )),
    print = FALSE
  )
  on.exit(untrace(EnhancedVolcano::EnhancedVolcano), add = TRUE)
  on.exit(options(mosuite_test_volcano_label_args = NULL), add = TRUE)

  boundary_data <- data.frame(
    Gene = c("both_boundary", "p_boundary", "fc_boundary", "neither"),
    annotation = c("a", "b", "c", "d"),
    `B-A_logFC` = c(1, 0.5, 1, 0.5),
    `B-A_pval` = c(0.05, 0.05, 0.06, 0.06),
    check.names = FALSE
  )

  plot_volcano_enhanced(
    boundary_data,
    feature_id_colname = "Gene",
    change_colname = "B-A_logFC",
    signif_colname = "B-A_pval",
    label_significant_features_only = TRUE,
    save_plots = FALSE,
    print_plots = FALSE
  )

  captured_labels <- getOption("mosuite_test_volcano_label_args")[[1]]$selectLab
  expect_length(captured_labels, 0)
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
  expect_equal(
    captured_args$col,
    c("grey30", "forestgreen", "royalblue", "red2")
  )
  expect_equal(captured_args$cutoffLineCol, "black")
})

test_that("plot_volcano_enhanced defaults match summary-style text sizing", {
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
          pointSize = pointSize
        ))
      )
    )),
    print = FALSE
  )
  on.exit(untrace(EnhancedVolcano::EnhancedVolcano), add = TRUE)
  on.exit(options(mosuite_test_volcano_args = NULL), add = TRUE)

  result <- plot_volcano_enhanced(
    nidap_deg_analysis,
    change_colname = "B-A_logFC",
    signif_colname = "B-A_pval",
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  captured_args <- getOption("mosuite_test_volcano_args")[[1]]
  expect_equal(captured_args$axisLabSize, 24)
  expect_equal(captured_args$titleLabSize, 24)
  expect_equal(captured_args$labSize, 7)
  expect_false(captured_args$drawConnectors)
  expect_equal(captured_args$pointSize, 2)

  volcano_plot <- attr(result, "plots")[[1]]
  expect_equal(volcano_plot$theme$axis.text$size, 16)
})

test_that("plot_volcano_enhanced supports custom axis labels and padding", {
  options(mosuite_test_volcano_args = list())
  trace(
    EnhancedVolcano::EnhancedVolcano,
    tracer = quote(options(
      mosuite_test_volcano_args = append(
        getOption("mosuite_test_volcano_args"),
        list(list(
          xlab = xlab,
          ylab = ylab,
          xlim = xlim,
          ylim = ylim
        ))
      )
    )),
    print = FALSE
  )
  on.exit(untrace(EnhancedVolcano::EnhancedVolcano), add = TRUE)
  on.exit(options(mosuite_test_volcano_args = NULL), add = TRUE)

  volcano_data <- data.frame(
    Gene = c("left_gene", "right_gene", "not_significant"),
    `B-A_logFC` = c(-1.2, 2.4, 0.2),
    `B-A_pval` = c(0.001, 0.02, 0.5),
    check.names = FALSE
  )

  result <- plot_volcano_enhanced(
    volcano_data,
    feature_id_colname = "Gene",
    change_colname = "B-A_logFC",
    signif_colname = "B-A_pval",
    change_lfc_name = "Custom fold change",
    change_sig_name = "Custom significance",
    use_custom_lab = TRUE,
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  captured_args <- getOption("mosuite_test_volcano_args")[[1]]
  expect_equal(captured_args$xlab, "Custom fold change")
  expect_equal(captured_args$ylab, "Custom significance")
  expect_equal(captured_args$xlim, c(-2, 3))
  expect_equal(captured_args$ylim, c(0, 3))
})

test_that("plot_volcano_enhanced matches selected-gene summary styling", {
  options(mosuite_test_volcano_args = list())
  trace(
    EnhancedVolcano::EnhancedVolcano,
    tracer = quote(options(
      mosuite_test_volcano_args = append(
        getOption("mosuite_test_volcano_args"),
        list(list(
          title = title,
          subtitle = subtitle,
          selectLab = selectLab,
          labSize = labSize,
          labCol = labCol
        ))
      )
    )),
    print = FALSE
  )
  on.exit(untrace(EnhancedVolcano::EnhancedVolcano), add = TRUE)
  on.exit(options(mosuite_test_volcano_args = NULL), add = TRUE)

  volcano_data <- data.frame(
    Gene = c("left_gene", "right_gene", "not_selected"),
    `G2-G1_logFC` = c(-1.4, 1.5, 0.2),
    `G2-G1_adjpval` = c(0.01, 0.02, 0.5),
    `G2-G1_tstat` = c(-10, 9, 1),
    check.names = FALSE
  )
  selected_genes <- c("left_gene", "right_gene")

  result <- plot_volcano_enhanced(
    volcano_data,
    feature_id_colname = "Gene",
    signif_colname = "G2-G1_adjpval",
    change_colname = "G2-G1_logFC",
    signif_threshold = 0.05,
    change_threshold = 1,
    value_to_sort_the_output_dataset = "t-statistic",
    num_features_to_label = 20,
    label_features = TRUE,
    custom_gene_list = paste(selected_genes, collapse = ","),
    custom_label_color = "black",
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  captured_args <- getOption("mosuite_test_volcano_args")[[1]]
  expect_equal(captured_args$title, "G2-G1")
  expect_null(captured_args$subtitle)
  expect_setequal(captured_args$selectLab, selected_genes)
  expect_equal(captured_args$labSize, 7)
  expect_true(all(captured_args$labCol == "black"))
})

test_that("plot_volcano_enhanced mimics summary contrast title and t-statistic labels", {
  options(mosuite_test_volcano_args = list())
  trace(
    EnhancedVolcano::EnhancedVolcano,
    tracer = quote(options(
      mosuite_test_volcano_args = append(
        getOption("mosuite_test_volcano_args"),
        list(list(
          title = title,
          subtitle = subtitle,
          selectLab = selectLab
        ))
      )
    )),
    print = FALSE
  )
  on.exit(untrace(EnhancedVolcano::EnhancedVolcano), add = TRUE)
  on.exit(options(mosuite_test_volcano_args = NULL), add = TRUE)

  volcano_data <- data.frame(
    Gene = c("p_value_top", "tstat_top", "fold_change_top"),
    `G2-G1_logFC` = c(1.2, -1.3, 2.5),
    `G2-G1_pval` = c(0.0001, 0.02, 0.03),
    `G2-G1_tstat` = c(2, -10, 5),
    check.names = FALSE
  )

  result <- plot_volcano_enhanced(
    volcano_data,
    feature_id_colname = "Gene",
    change_colname = "G2-G1_logFC",
    signif_colname = "G2-G1_pval",
    value_to_sort_the_output_dataset = "t-statistic",
    num_features_to_label = 1,
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(result, "data.frame")
  captured_args <- getOption("mosuite_test_volcano_args")[[1]]
  expect_equal(captured_args$title, "G2-G1")
  expect_null(captured_args$subtitle)
  expect_equal(captured_args$selectLab, "tstat_top")
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

test_that("plot_volcano_enhanced offsets labels when connectors are enabled", {
  volcano_data <- data.frame(
    Gene = c("gene_a", "gene_b", "not_selected"),
    `B-A_logFC` = c(2, -2, 0.1),
    `B-A_pval` = c(0.001, 0.002, 0.5),
    check.names = FALSE
  )

  result <- suppressWarnings(plot_volcano_enhanced(
    volcano_data,
    feature_id_colname = "Gene",
    change_colname = "B-A_logFC",
    signif_colname = "B-A_pval",
    label_features = TRUE,
    custom_gene_list = "gene_a,gene_b",
    draw_connectors = TRUE,
    save_plots = FALSE,
    print_plots = FALSE
  ))

  expect_s3_class(result, "data.frame")
  volcano_plot <- attr(result, "plots")[[1]]
  layer_geoms <- vapply(
    volcano_plot$layers,
    function(layer) class(layer$geom)[1],
    character(1)
  )
  expect_true("GeomTextRepel" %in% layer_geoms)
  expect_false("GeomText" %in% layer_geoms)
})

test_that("plot_volcano_enhanced preserves filename for one comparison", {
  options(mosuite_test_plot_output_args = list())
  trace(
    ggplot2::ggsave,
    tracer = quote({
      options(
        mosuite_test_plot_output_args = append(
          getOption("mosuite_test_plot_output_args"),
          list(list(filename = filename))
        )
      )
    }),
    print = FALSE
  )
  on.exit(untrace(ggplot2::ggsave), add = TRUE)
  on.exit(options(mosuite_test_plot_output_args = NULL), add = TRUE)

  plots_dir <- tempfile("volcano-single-output-")
  dir.create(plots_dir)
  on.exit(unlink(plots_dir, recursive = TRUE), add = TRUE)

  expect_no_error(
    result <- plot_volcano_enhanced(
      nidap_deg_analysis,
      change_colname = "B-A_logFC",
      signif_colname = "B-A_pval",
      save_plots = TRUE,
      print_plots = FALSE,
      plots_subdir = plots_dir,
      plot_filename = "custom_volcano.png"
    )
  )

  expect_s3_class(result, "data.frame")
  captured_output_args <- getOption("mosuite_test_plot_output_args")
  expect_length(captured_output_args, 1)
  expect_match(
    captured_output_args[[1]]$filename,
    basename(plots_dir),
    fixed = TRUE
  )
  expect_match(
    captured_output_args[[1]]$filename,
    "custom_volcano.png",
    fixed = TRUE
  )
})

test_that("plot_volcano_enhanced auto-detects change_colname and signif_colname when NULL", {
  # Regression test: all three contrasts in nidap_deg_analysis should be
  # detected automatically, preferring _adjpval over _pval.
  expect_no_error(
    result <- plot_volcano_enhanced(
      nidap_deg_analysis,
      save_plots = FALSE,
      print_plots = FALSE
    )
  )
  expect_s3_class(result, "data.frame")
  expect_length(
    attr(result, "plots"),
    length(grep("_logFC$", colnames(nidap_deg_analysis)))
  )
})

test_that("plot_volcano_enhanced errors when no _logFC columns and change_colname is NULL", {
  bad_df <- data.frame(
    Gene = letters[1:3],
    `B-A_adjpval` = c(0.01, 0.02, 0.5),
    check.names = FALSE
  )
  expect_error(
    plot_volcano_enhanced(bad_df, save_plots = FALSE, print_plots = FALSE),
    regexp = "_logFC"
  )
})

test_that("plot_volcano_enhanced errors when no adjpval/pval columns and signif_colname is NULL", {
  bad_df <- data.frame(
    Gene = letters[1:3],
    `B-A_logFC` = c(1, -1, 0.1),
    check.names = FALSE
  )
  expect_error(
    plot_volcano_enhanced(bad_df, save_plots = FALSE, print_plots = FALSE),
    regexp = "auto-detect"
  )
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
