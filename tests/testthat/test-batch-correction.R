test_that("batch_correction works for NIDAP", {
  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts),
      "norm" = list("voom" = as.data.frame(nidap_norm_counts))
    )
  ) |>
    batch_correct_counts(
      count_type = "norm",
      sub_count_type = "voom",
      covariates_colnames = "Group",
      batch_colname = "Batch",
      label_colname = "Label",
      print_plots = TRUE
    )
  # TODO: getting different results than nidap_batch_corrected_counts
  expect_true(all.equal(
    moo@counts[["batch"]] |>
      dplyr::arrange(desc(Gene)),
    as.data.frame(nidap_batch_corrected_counts_2) |>
      dplyr::arrange(desc(Gene))
  ))
})

test_that("batch_correction warnings & errors", {
  moo <- create_multiOmicDataSet_from_dataframes(
    readr::read_tsv(
      system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite")
    ) |>
      dplyr::mutate(batch = 1),
    gene_counts
  ) |>
    clean_raw_counts() |>
    filter_counts(
      group_colname = "condition",
      label_colname = "sample_id",
      minimum_count_value_to_be_considered_nonzero = 1,
      minimum_number_of_samples_with_nonzero_counts_in_total = 1,
      minimum_number_of_samples_with_nonzero_counts_in_a_group = 1,
      print_plots = FALSE
    ) |>
    normalize_counts(group_colname = "condition", label_colname = "sample_id")

  expect_warning(
    moo |>
      batch_correct_counts(
        covariates_colnames = "condition",
        batch_colname = "batch",
        label_colname = "sample_id"
      ),
    "Batch column 'batch' contains only 1 unique value"
  )
  expect_error(
    moo |>
      batch_correct_counts(
        covariates_colnames = "batch",
        batch_colname = "batch"
      ),
    "Batch column 'batch' cannot be included in covariates."
  )
})

test_that("batch_correct_counts forwards plot settings to PCA and histogram", {
  pca_args <- NULL
  histogram_args <- NULL

  local_mocked_bindings(
    ComBat = function(dat, ...) dat,
    .package = "sva"
  )
  local_mocked_bindings(
    plot_pca = function(...) {
      pca_args <<- list(...)
      return(ggplot2::ggplot())
    },
    plot_histogram = function(...) {
      histogram_args <<- list(...)
      return(ggplot2::ggplot())
    },
    print_or_save_plot = function(...) invisible(NULL),
    .package = "MOSuite"
  )

  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts),
      "norm" = list("voom" = as.data.frame(nidap_norm_counts))
    )
  )

  batch_correct_counts(
    moo,
    count_type = "norm",
    sub_count_type = "voom",
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    covariates_colnames = "Group",
    batch_colname = "Batch",
    label_colname = "Label",
    samples_to_rename = c("A1:Alpha 1"),
    add_label_to_pca = FALSE,
    principal_component_on_x_axis = 2,
    principal_component_on_y_axis = 3,
    legend_position_for_pca = "bottom",
    label_offset_x_ = 4,
    label_offset_y_ = 5,
    label_font_size = 6,
    point_size_for_pca = 7,
    color_histogram_by_group = FALSE,
    set_min_max_for_x_axis_for_histogram = TRUE,
    minimum_for_x_axis_for_histogram = -2,
    maximum_for_x_axis_for_histogram = 2,
    legend_font_size_for_histogram = 11,
    legend_position_for_histogram = "right",
    number_of_histogram_legend_columns = 2,
    colors_for_plots = c(A = "red", B = "blue", C = "green"),
    plot_corr_matrix_heatmap = FALSE,
    print_plots = TRUE,
    save_plots = FALSE
  )

  expect_equal(pca_args$samples_to_rename, c("A1:Alpha 1"))
  expect_equal(pca_args$principal_components, c(2, 3))
  expect_equal(pca_args$legend_position, "bottom")
  expect_equal(pca_args$point_size, 7)
  expect_null(pca_args$label_colname)
  expect_equal(pca_args$label_font_size, 6)
  expect_equal(pca_args$label_offset_x_, 4)
  expect_equal(pca_args$label_offset_y_, 5)
  expect_equal(pca_args$color_values, c(A = "red", B = "blue", C = "green"))

  expect_false(histogram_args$color_by_group)
  expect_true(histogram_args$set_min_max_for_x_axis)
  expect_equal(histogram_args$minimum_for_x_axis, -2)
  expect_equal(histogram_args$maximum_for_x_axis, 2)
  expect_equal(histogram_args$x_axis_label, "Batch Corrected Counts")
  expect_equal(histogram_args$legend_font_size, 11)
  expect_equal(histogram_args$legend_position, "right")
  expect_equal(histogram_args$number_of_legend_columns, 2)
  expect_equal(histogram_args$color_values, moo@analyses[["colors"]][["Label"]])

  pca_args <- NULL
  batch_correct_counts(
    moo,
    count_type = "norm",
    sub_count_type = "voom",
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    covariates_colnames = "Group",
    batch_colname = "Batch",
    label_colname = "Label",
    add_label_to_pca = TRUE,
    plot_corr_matrix_heatmap = FALSE,
    print_plots = TRUE,
    save_plots = FALSE
  )
  expect_equal(pca_args$label_colname, "Label")
})

test_that("batch_correct_counts handles histogram label combinations", {
  pca_args <- NULL
  histogram_args <- NULL
  group_colors <- c(A = "red", B = "blue", C = "green")

  local_mocked_bindings(
    ComBat = function(dat, ...) dat,
    .package = "sva"
  )
  local_mocked_bindings(
    plot_pca = function(...) {
      pca_args <<- list(...)
      return(ggplot2::ggplot())
    },
    plot_histogram = function(...) {
      histogram_args <<- list(...)
      return(ggplot2::ggplot())
    },
    print_or_save_plot = function(...) invisible(NULL),
    .package = "MOSuite"
  )

  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts),
      "norm" = list("voom" = as.data.frame(nidap_norm_counts))
    )
  )

  combinations <- list(
    list(
      label_colname = NULL,
      color_histogram_by_group = FALSE,
      interactive_plots = FALSE
    ),
    list(
      label_colname = NULL,
      color_histogram_by_group = FALSE,
      interactive_plots = TRUE
    ),
    list(
      label_colname = NULL,
      color_histogram_by_group = TRUE,
      interactive_plots = FALSE
    ),
    list(
      label_colname = NULL,
      color_histogram_by_group = TRUE,
      interactive_plots = TRUE
    ),
    list(
      label_colname = "Label",
      color_histogram_by_group = FALSE,
      interactive_plots = FALSE
    ),
    list(
      label_colname = "Label",
      color_histogram_by_group = FALSE,
      interactive_plots = TRUE
    ),
    list(
      label_colname = "Label",
      color_histogram_by_group = TRUE,
      interactive_plots = FALSE
    ),
    list(
      label_colname = "Label",
      color_histogram_by_group = TRUE,
      interactive_plots = TRUE
    )
  )

  for (combination in combinations) {
    pca_args <- NULL
    histogram_args <- NULL
    batch_correct_counts(
      moo,
      count_type = "norm",
      sub_count_type = "voom",
      sample_id_colname = "Sample",
      feature_id_colname = "Gene",
      covariates_colnames = "Group",
      batch_colname = "Batch",
      label_colname = combination$label_colname,
      color_histogram_by_group = combination$color_histogram_by_group,
      interactive_plots = combination$interactive_plots,
      colors_for_plots = group_colors,
      plot_corr_matrix_heatmap = FALSE,
      print_plots = TRUE,
      save_plots = FALSE
    )

    expected_label_colname <- if (is.null(combination$label_colname)) {
      "Sample"
    } else {
      combination$label_colname
    }
    expected_histogram_colors <- if (
      isTRUE(combination$color_histogram_by_group)
    ) {
      group_colors
    } else {
      moo@analyses[["colors"]][[expected_label_colname]]
    }

    expect_equal(pca_args$label_colname, combination$label_colname)
    expect_equal(histogram_args$label_colname, expected_label_colname)
    expect_equal(
      histogram_args$color_by_group,
      combination$color_histogram_by_group
    )
    expect_equal(
      histogram_args$interactive_plots,
      combination$interactive_plots
    )
    expect_equal(histogram_args$color_values, expected_histogram_colors)
  }
})

test_that("batch_correct_counts forwards the default MOSuite plot colors", {
  pca_args <- NULL
  histogram_args <- NULL
  expected_colors <- c(
    "1" = "#ff9287",
    "2" = "#008cf9"
  )

  local_mocked_bindings(
    ComBat = function(dat, ...) dat,
    .package = "sva"
  )
  local_mocked_bindings(
    plot_pca = function(...) {
      pca_args <<- list(...)
      return(ggplot2::ggplot())
    },
    plot_histogram = function(...) {
      histogram_args <<- list(...)
      return(ggplot2::ggplot())
    },
    print_or_save_plot = function(...) invisible(NULL),
    .package = "MOSuite"
  )

  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts),
      "norm" = list("voom" = as.data.frame(nidap_norm_counts))
    )
  )

  batch_correct_counts(
    moo,
    count_type = "norm",
    sub_count_type = "voom",
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    covariates_colnames = "Group",
    batch_colname = "Batch",
    label_colname = "Label",
    plot_corr_matrix_heatmap = FALSE,
    print_plots = TRUE,
    save_plots = FALSE
  )

  expect_equal(pca_args$color_values, expected_colors)
  expect_equal(histogram_args$color_values, expected_colors)
})

test_that("batch_correct_counts PCA matches standalone plot_pca on batch output", {
  pca_capture <- capture_saved_pca_plot()
  local_mocked_bindings(
    print_or_save_plot = pca_capture$print_or_save_plot,
    .package = "MOSuite"
  )

  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts),
      "norm" = list("voom" = as.data.frame(nidap_norm_counts))
    )
  ) |>
    batch_correct_counts(
      count_type = "norm",
      sub_count_type = "voom",
      sample_id_colname = "Sample",
      feature_id_colname = "Gene",
      covariates_colnames = "Group",
      batch_colname = "Batch",
      label_colname = NULL,
      plot_corr_matrix_heatmap = FALSE,
      print_plots = TRUE,
      save_plots = FALSE
    )

  expected_pca <- plot_pca(
    moo@counts$batch,
    sample_metadata = moo@sample_meta,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    group_colname = "Batch",
    label_colname = NULL,
    samples_to_rename = c(""),
    principal_components = c(1, 2),
    legend_position = "top",
    point_size = 5,
    label_font_size = 3,
    label_offset_y_ = 2,
    label_offset_x_ = 2,
    log_transform = FALSE,
    print_plots = FALSE,
    save_plots = FALSE
  )

  expect_s3_class(pca_capture$get(), "ggplot")
  expect_pca_coordinates_equal(pca_capture$get(), expected_pca)
})

test_that("batch_correct_counts histogram matches standalone plot_histogram on batch output", {
  histogram_capture <- capture_saved_histogram_plot()
  local_mocked_bindings(
    print_or_save_plot = histogram_capture$print_or_save_plot,
    .package = "MOSuite"
  )

  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts),
      "norm" = list("voom" = as.data.frame(nidap_norm_counts))
    )
  ) |>
    batch_correct_counts(
      count_type = "norm",
      sub_count_type = "voom",
      sample_id_colname = "Sample",
      feature_id_colname = "Gene",
      covariates_colnames = "Group",
      batch_colname = "Batch",
      label_colname = NULL,
      plot_corr_matrix_heatmap = FALSE,
      print_plots = TRUE,
      save_plots = FALSE,
      interactive_plots = FALSE
    )

  expected_histogram <- plot_histogram(
    moo@counts$batch,
    sample_metadata = moo@sample_meta,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    group_colname = "Batch",
    color_values = moo@analyses$colors[["Batch"]],
    label_colname = NULL,
    color_by_group = TRUE,
    set_min_max_for_x_axis = FALSE,
    minimum_for_x_axis = -1,
    maximum_for_x_axis = 1,
    x_axis_label = "Batch Corrected Counts",
    legend_position = "top",
    legend_font_size = NULL,
    number_of_legend_columns = 6,
    interactive_plots = FALSE
  ) +
    ggplot2::labs(caption = "batch-corrected counts")

  expect_s3_class(histogram_capture$get(), "ggplot")
  expect_histogram_layers_equal(histogram_capture$get(), expected_histogram)
})
