test_that("normalize works for NIDAP", {
  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts)
    )
  ) |>
    normalize_counts(
      group_colname = "Group",
      label_colname = "Label",
      print_plots = TRUE
    )
  expect_true(equal_dfs(
    moo@counts[["norm"]][["voom"]] |>
      dplyr::arrange(desc(Gene)),
    as.data.frame(nidap_norm_counts) |>
      dplyr::arrange(desc(Gene))
  ))
})

test_that("normalize works for RENEE", {
  moo <- create_multiOmicDataSet_from_dataframes(
    readr::read_tsv(
      system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite")
    ),
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
  expect_equal(
    head(moo@counts$norm$voom),
    structure(
      list(
        gene_id = c(
          "ENSG00000215458.8",
          "ENSG00000160179.18",
          "ENSG00000258017.1",
          "ENSG00000282393.1",
          "ENSG00000286104.1",
          "ENSG00000274422.1"
        ),
        KO_S3 = c(
          11.0751960068561,
          9.6086338540783,
          9.6086338540783,
          8.81615260371772,
          9.6086338540783,
          8.81615260371772
        ),
        KO_S4 = c(
          12.3480907442867,
          12.7703165561761,
          8.81615260371772,
          9.6086338540783,
          8.81615260371772,
          9.6086338540783
        ),
        WT_S1 = c(
          8.81615260371772,
          12.3480907442867,
          8.81615260371772,
          8.81615260371772,
          8.81615260371772,
          8.81615260371772
        ),
        WT_S2 = c(
          10.0048744792586,
          12.2369960496953,
          8.81615260371772,
          8.81615260371772,
          8.81615260371772,
          8.81615260371772
        )
      ),
      row.names = c(NA, 6L),
      class = "data.frame"
    )
  )
  expect_equal(
    tail(moo@counts$norm$voom),
    structure(
      list(
        gene_id = c(
          "ENSG00000157538.14",
          "ENSG00000160193.11",
          "ENSG00000182093.15",
          "ENSG00000182362.14",
          "ENSG00000173276.14",
          "ENSG00000237232.7"
        ),
        KO_S3 = c(
          12.3480907442867,
          9.6086338540783,
          11.8597009422769,
          11.0751960068561,
          11.8597009422769,
          8.81615260371772
        ),
        KO_S4 = c(
          12.7703165561761,
          9.6086338540783,
          9.6086338540783,
          8.81615260371772,
          12.7703165561761,
          9.6086338540783
        ),
        WT_S1 = c(
          12.2426956580003,
          10.5853565029804,
          11.7865266999202,
          8.81615260371772,
          11.7865266999202,
          8.81615260371772
        ),
        WT_S2 = c(
          12.4720029602325,
          10.9249210977479,
          11.4357186116131,
          8.81615260371772,
          12.2369960496953,
          8.81615260371772
        )
      ),
      row.names = 286:291,
      class = "data.frame"
    )
  )
})

test_that("normalize_counts forwards plotting parameters", {
  pca_args <- NULL
  histogram_args <- NULL

  local_mocked_bindings(
    plot_pca = function(...) {
      pca_args <<- list(...)
      ggplot2::ggplot()
    },
    plot_histogram = function(...) {
      histogram_args <<- list(...)
      ggplot2::ggplot()
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
      "filt" = as.data.frame(nidap_filtered_counts)
    )
  )

  normalize_counts(
    moo,
    group_colname = "Group",
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
  expect_false(pca_args$add_label)
  expect_equal(pca_args$label_font_size, 6)
  expect_equal(pca_args$label_offset_x_, 4)
  expect_equal(pca_args$label_offset_y_, 5)
  expect_equal(pca_args$color_values, c(A = "red", B = "blue", C = "green"))

  expect_false(histogram_args$color_by_group)
  expect_true(histogram_args$set_min_max_for_x_axis)
  expect_equal(histogram_args$minimum_for_x_axis, -2)
  expect_equal(histogram_args$maximum_for_x_axis, 2)
  expect_equal(histogram_args$legend_font_size, 11)
  expect_equal(histogram_args$legend_position, "right")
  expect_equal(histogram_args$number_of_legend_columns, 2)
  expect_equal(histogram_args$color_values, moo@analyses[["colors"]][["Label"]])
})

test_that("normalize_counts forwards the default MOSuite plot colors", {
  pca_args <- NULL
  histogram_args <- NULL
  default_colors <- c(
    "#5954d6",
    "#e1562c",
    "#b80058",
    "#00c6f8",
    "#d163e6",
    "#00a76c",
    "#ff9287",
    "#008cf9",
    "#006e00",
    "#796880",
    "#FFA500",
    "#878500"
  )

  local_mocked_bindings(
    plot_pca = function(...) {
      pca_args <<- list(...)
      ggplot2::ggplot()
    },
    plot_histogram = function(...) {
      histogram_args <<- list(...)
      ggplot2::ggplot()
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
      "filt" = as.data.frame(nidap_filtered_counts)
    )
  )

  normalize_counts(
    moo,
    group_colname = "Group",
    label_colname = "Label",
    plot_corr_matrix_heatmap = FALSE,
    print_plots = TRUE,
    save_plots = FALSE
  )

  expect_equal(pca_args$color_values, default_colors)
  expect_equal(histogram_args$color_values, default_colors)
})
