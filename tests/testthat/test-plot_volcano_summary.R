test_that("plot_volcano_summary works on nidap dataset", {
  expect_snapshot(
    df_volc_sum <- plot_volcano_summary(
      nidap_deg_analysis,
      save_plots = FALSE,
      print_plots = FALSE
    )
  )
  expect_equal(
    head(df_volc_sum),
    structure(
      list(
        Gene = c("Dntt", "Tmsb4x", "Flt3", "Tspan13", "Tapt1", "Itgb7"),
        Contrast = c("B-A", "B-A", "B-A", "B-A", "B-A", "B-A"),
        FC = c(
          -42.7465863415622,
          3.85002020608143,
          -7.71439441748029,
          -7.03849783123801,
          -5.29181569343323,
          8.87382341151917
        ),
        logFC = c(
          -5.41773730869316,
          1.94486601753143,
          -2.94755290920186,
          -2.81526755916543,
          -2.40376281543362,
          3.14955584391085
        ),
        tstat = c(
          -15.6879749543426,
          12.9102607749226,
          -11.3808403447749,
          -11.0312744854072,
          -10.6584674633331,
          10.5614738819538
        ),
        pval = c(
          3.15934346857821e-09,
          2.76055502226637e-08,
          1.09340538530663e-07,
          1.53110956271563e-07,
          2.21459280934843e-07,
          2.44206995658642e-07
        ),
        adjpval = c(
          2.50946651709167e-05,
          0.000109635442709309,
          0.000289497299183018,
          0.000304040081416256,
          0.000323289361086099,
          0.000323289361086099
        )
      ),
      row.names = c("B-A.1", "B-A.2", "B-A.3", "B-A.4", "B-A.5", "B-A.6"),
      class = "data.frame"
    )
  )
  expect_equal(
    tail(df_volc_sum),
    structure(
      list(
        Gene = c("Tecpr1", "Lap3", "Zfp952", "Tsr3", "Nbas", "Slc50a1"),
        Contrast = c("B-C", "B-C", "B-C", "B-C", "B-C", "B-C"),
        FC = c(
          -17.6925615963148,
          2.57712293045075,
          -10.2589472087027,
          -3.22520189762021,
          4.43444692871868,
          -2.36807519790042
        ),
        logFC = c(
          -4.14507103690983,
          1.36576135633674,
          -3.35881078151314,
          -1.68938947606405,
          2.1487541805238,
          -1.24371489427577
        ),
        tstat = c(
          -2.19458425130448,
          2.1944817392618,
          -2.19287280238278,
          -2.19094226223025,
          2.19079653013039,
          -2.18921321212647
        ),
        pval = c(
          0.0491166107800282,
          0.0491255772255026,
          0.0492665099590459,
          0.0494361189306691,
          0.0494489447669875,
          0.0495884954830029
        ),
        adjpval = c(
          0.265806852794392,
          0.265806852794392,
          0.266387943229885,
          0.26682946214958,
          0.26682946214958,
          0.26740082798472
        )
      ),
      row.names = c(
        "B-C.957",
        "B-C.958",
        "B-C.959",
        "B-C.960",
        "B-C.961",
        "B-C.962"
      ),
      class = "data.frame"
    )
  )
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

  result <- plot_volcano_summary(
    nidap_deg_analysis,
    add_features = TRUE,
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
