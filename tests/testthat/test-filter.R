test_that("filter_counts reproduces NIDAP results", {
  set.seed(10)
  moo <- create_multiOmicDataSet_from_dataframes(
    as.data.frame(nidap_sample_metadata),
    as.data.frame(nidap_clean_raw_counts),
    sample_id_colname = "Sample",
    feature_id_colname = "Gene"
  ) |>
    calc_cpm(feature_id_colname = "Gene") |>
    filter_counts(
      sample_id_colname = "Sample",
      feature_id_colname = "Gene",
      count_type = "raw",
      print_plots = TRUE
    )
  rds_counts_filt <- moo@counts$filt |>
    dplyr::arrange(desc(Gene))
  nidap_counts_filt <- as.data.frame(nidap_filtered_counts) |>
    dplyr::arrange(desc(Gene))

  expect_true(equal_dfs(rds_counts_filt, nidap_counts_filt))
})

# TODO get filter_counts() to work on tibbles too, not only dataframes

test_that("filter_counts works on RENEE dataset", {
  moo <- create_multiOmicDataSet_from_dataframes(
    readr::read_tsv(
      system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite")
    ),
    gene_counts |> glue_gene_symbols()
  )
  rds2 <- moo |>
    filter_counts(
      feature_id_colname = "gene_id",
      sample_id_colname = "sample_id",
      group_colname = "condition",
      label_colname = "sample_id",
      samples_to_include = c("KO_S3", "KO_S4", "WT_S1", "WT_S2"),
      minimum_count_value_to_be_considered_nonzero = 1,
      minimum_number_of_samples_with_nonzero_counts_in_total = 1,
      minimum_number_of_samples_with_nonzero_counts_in_a_group = 1,
      print_plots = TRUE,
      count_type = "raw"
    )
  expect_equal(dim(rds2@counts$filt), c(291, 5))
  expect_equal(
    rds2@counts$filt |> dplyr::arrange(gene_id) |> head(),
    structure(
      list(
        gene_id = c(
          "ENSG00000072803.17|FBXW11",
          "ENSG00000083845.9|RPS5",
          "ENSG00000107371.13|EXOSC3",
          "ENSG00000111639.8|MRPL51",
          "ENSG00000111640.15|GAPDH",
          "ENSG00000111786.9|SRSF9"
        ),
        KO_S3 = c(2, 1, 1, 0, 0, 0),
        KO_S4 = c(0, 0, 1, 1, 1, 1),
        WT_S1 = c(0, 0, 0, 0, 0, 0),
        WT_S2 = c(0, 0, 0, 0, 0, 0)
      ),
      row.names = c(NA, 6L),
      class = "data.frame"
    )
  )
  expect_equal(
    rds2@counts$filt |> dplyr::arrange(gene_id) |> tail(),
    structure(
      list(
        gene_id = c(
          "ENSG00000281903.2|LINC02246",
          "ENSG00000282393.1|AC016588.2",
          "ENSG00000283886.2|BX664615.2",
          "ENSG00000285413.1|AP001056.2",
          "ENSG00000286018.1|AF129075.3",
          "ENSG00000286104.1|AC016629.3"
        ),
        KO_S3 = c(0.85, 0, 1, 3, 2, 1),
        KO_S4 = c(0, 1, 0, 1, 0, 0),
        WT_S1 = c(0, 0, 0, 0, 0, 0),
        WT_S2 = c(0.71, 0, 0, 0, 0, 0)
      ),
      row.names = 286:291,
      class = "data.frame"
    )
  )
})

test_that("remove_low_count_genes works", {
  df <- data.frame(
    Gene = c(
      "mt-Nd5_43275",
      "mt-Nd6_43276",
      "mt-Te_43277",
      "mt-Cytb_43278",
      "mt-Tt_43279",
      "mt-Tp_43280"
    ),
    A1 = c(6155, 858, 0, 20542, 0, 12),
    A2 = c(10823, 1420, 1, 29677, 9, 16),
    A3 = c(9482, 1167, 2, 31730, 0, 13),
    B1 = c(6162, 1181, 0, 28293, 0, 15),
    B2 = c(8002, 845, 1, 25617, 7, 19),
    B3 = c(7225, 1198, 3, 30370, 3, 26),
    C1 = c(4141, 515, 0, 21310, 0, 32),
    C2 = c(9058, 1147, 4, 30108, 0, 33),
    C3 = c(8481, 1124, 2, 30893, 2, 50),
    row.names = seq(5, 10)
  )
  sample_meta <- structure(
    list(
      Sample = c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"),
      Group = c("A", "A", "A", "B", "B", "B", "C", "C", "C"),
      Replicate = c(1, 2, 3, 1, 2, 3, 1, 2, 3),
      Batch = c(1, 2, 2, 1, 1, 2, 1, 2, 2),
      Label = c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")
    ),
    row.names = c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"),
    class = "data.frame"
  )

  # test default params
  expect_equal(
    remove_low_count_genes(
      counts_dat = df,
      sample_metadata = sample_meta,
      feature_id_colname = "Gene",
      group_colname = "Group",
      use_cpm_counts_to_filter = TRUE,
      use_group_based_filtering = FALSE,
      minimum_count_value_to_be_considered_nonzero = 8,
      minimum_number_of_samples_with_nonzero_counts_in_total = 7,
      minimum_number_of_samples_with_nonzero_counts_in_a_group = 3
    ),
    structure(
      list(
        Gene = c(
          "mt-Nd5_43275",
          "mt-Nd6_43276",
          "mt-Cytb_43278",
          "mt-Tp_43280"
        ),
        A1 = c(
          223274.204664998,
          31124.1702035042,
          745166.322051728,
          435.303079769289
        ),
        A2 = c(
          258022.219043532,
          33853.0491584418,
          707504.88723597,
          381.442807419063
        ),
        A3 = c(
          223663.725998962,
          27527.4803038166,
          748454.970042931,
          306.647167051941
        ),
        B1 = c(
          172842.276513983,
          33126.7005133096,
          793610.277411573,
          420.74556113433
        ),
        B2 = c(
          232002.551390218,
          24499.1447044156,
          742715.490997652,
          550.868342466151
        ),
        B3 = c(
          186091.435930457,
          30856.406954282,
          782227.94591114,
          669.671603348358
        ),
        C1 = c(
          159281.483191015,
          19809.2160935457,
          819678.436802831,
          1230.86391260866
        ),
        C2 = c(
          224485.749690211,
          28426.2701363073,
          746171.003717472,
          817.843866171004
        ),
        C3 = c(
          209138.883408956,
          27717.4985204182,
          761811.994476228,
          1232.98480962715
        )
      ),
      row.names = 1:4,
      class = "data.frame"
    )
  )
})

test_that("remove_low_count_genes counts threshold values inclusively", {
  df <- data.frame(
    Gene = c("exact_keep", "below_minimum", "above_keep"),
    S1 = c(4, 4, 5),
    S2 = c(5, 0, 5),
    S3 = c(0, 0, 0),
    check.names = FALSE
  )
  sample_meta <- data.frame(
    Sample = c("S1", "S2", "S3"),
    Group = c("A", "A", "A"),
    row.names = c("S1", "S2", "S3"),
    check.names = FALSE
  )

  result <- remove_low_count_genes(
    counts_dat = df,
    sample_metadata = sample_meta,
    feature_id_colname = "Gene",
    group_colname = "Group",
    use_cpm_counts_to_filter = FALSE,
    use_group_based_filtering = FALSE,
    minimum_count_value_to_be_considered_nonzero = 4,
    minimum_number_of_samples_with_nonzero_counts_in_total = 2,
    minimum_number_of_samples_with_nonzero_counts_in_a_group = 2
  )

  expect_equal(result$Gene, c("exact_keep", "above_keep"))
})

test_that("remove_low_count_genes uses original filter threshold parameters", {
  df <- data.frame(
    Gene = c("keep", "remove"),
    S1 = c(4, 4),
    S2 = c(4, 0),
    check.names = FALSE
  )
  sample_meta <- data.frame(
    Sample = c("S1", "S2"),
    Group = c("A", "A"),
    row.names = c("S1", "S2"),
    check.names = FALSE
  )

  result <- remove_low_count_genes(
    counts_dat = df,
    sample_metadata = sample_meta,
    feature_id_colname = "Gene",
    group_colname = "Group",
    use_cpm_counts_to_filter = FALSE,
    use_group_based_filtering = FALSE,
    minimum_count_value_to_be_considered_nonzero = 4,
    minimum_number_of_samples_with_nonzero_counts_in_total = 2,
    minimum_number_of_samples_with_nonzero_counts_in_a_group = 2
  )

  expect_equal(result$Gene, "keep")
})

test_that("remove_low_count_genes ignores non-numeric annotation columns", {
  df <- data.frame(
    Gene = c("keep", "remove"),
    Symbol = c("A", "B"),
    S1 = c(4, 4),
    S2 = c(4, 0),
    check.names = FALSE
  )
  sample_meta <- data.frame(
    Sample = c("S1", "S2"),
    Group = c("A", "A"),
    row.names = c("S1", "S2"),
    check.names = FALSE
  )

  result <- remove_low_count_genes(
    counts_dat = df,
    sample_metadata = sample_meta,
    feature_id_colname = "Gene",
    group_colname = "Group",
    use_cpm_counts_to_filter = FALSE,
    use_group_based_filtering = FALSE,
    minimum_count_value_to_be_considered_nonzero = 4,
    minimum_number_of_samples_with_nonzero_counts_in_total = 2,
    minimum_number_of_samples_with_nonzero_counts_in_a_group = 2
  )

  expect_equal(result$Gene, "keep")
  expect_equal(result$Symbol, "A")
})

test_that("remove_low_count_genes distinguishes total sample minimum 8 from 9", {
  sample_names <- paste0("S", 1:10)
  df <- data.frame(
    Gene = c("passes_exactly_eight", "passes_nine"),
    S1 = c(4, 4),
    S2 = c(4, 4),
    S3 = c(4, 4),
    S4 = c(4, 4),
    S5 = c(4, 4),
    S6 = c(4, 4),
    S7 = c(4, 4),
    S8 = c(4, 4),
    S9 = c(0, 4),
    S10 = c(0, 0),
    check.names = FALSE
  )
  sample_meta <- data.frame(
    Sample = sample_names,
    Group = rep("A", length(sample_names)),
    row.names = sample_names,
    check.names = FALSE
  )

  total_8 <- remove_low_count_genes(
    counts_dat = df,
    sample_metadata = sample_meta,
    feature_id_colname = "Gene",
    group_colname = "Group",
    use_cpm_counts_to_filter = FALSE,
    use_group_based_filtering = FALSE,
    minimum_count_value_to_be_considered_nonzero = 4,
    minimum_number_of_samples_with_nonzero_counts_in_total = 8,
    minimum_number_of_samples_with_nonzero_counts_in_a_group = 2
  )
  total_9 <- remove_low_count_genes(
    counts_dat = df,
    sample_metadata = sample_meta,
    feature_id_colname = "Gene",
    group_colname = "Group",
    use_cpm_counts_to_filter = FALSE,
    use_group_based_filtering = FALSE,
    minimum_count_value_to_be_considered_nonzero = 4,
    minimum_number_of_samples_with_nonzero_counts_in_total = 9,
    minimum_number_of_samples_with_nonzero_counts_in_a_group = 2
  )

  expect_equal(total_8$Gene, c("passes_exactly_eight", "passes_nine"))
  expect_equal(total_9$Gene, "passes_nine")
})

test_that("remove_low_count_genes applies sample minimum within groups", {
  df <- data.frame(
    Gene = c(
      "group_exact_keep",
      "split_across_groups",
      "one_group_keep",
      "single_sample_total_only"
    ),
    S1 = c(4, 4, 0, 5),
    S2 = c(4, 0, 0, 0),
    S3 = c(0, 4, 5, 0),
    S4 = c(0, 0, 4, 0),
    check.names = FALSE
  )
  sample_meta <- data.frame(
    Sample = c("S1", "S2", "S3", "S4"),
    Group = c("A", "A", "B", "B"),
    row.names = paste0("metadata_row_", 1:4),
    check.names = FALSE
  )

  group_result <- remove_low_count_genes(
    counts_dat = df,
    sample_metadata = sample_meta,
    feature_id_colname = "Gene",
    group_colname = "Group",
    use_cpm_counts_to_filter = FALSE,
    use_group_based_filtering = TRUE,
    minimum_count_value_to_be_considered_nonzero = 4,
    minimum_number_of_samples_with_nonzero_counts_in_total = 1,
    minimum_number_of_samples_with_nonzero_counts_in_a_group = 2
  )
  non_group_result <- remove_low_count_genes(
    counts_dat = df,
    sample_metadata = sample_meta,
    feature_id_colname = "Gene",
    group_colname = "Group",
    use_cpm_counts_to_filter = FALSE,
    use_group_based_filtering = FALSE,
    minimum_count_value_to_be_considered_nonzero = 4,
    minimum_number_of_samples_with_nonzero_counts_in_total = 2,
    minimum_number_of_samples_with_nonzero_counts_in_a_group = 2
  )

  expect_equal(group_result$Gene, c("group_exact_keep", "one_group_keep"))
  expect_equal(
    non_group_result$Gene,
    c("group_exact_keep", "split_across_groups", "one_group_keep")
  )
})

test_that("remove_low_count_genes ignores missing group assignments", {
  df <- data.frame(
    Gene = c("keep", "remove"),
    S1 = c(4, 4),
    S2 = c(4, 0),
    S3 = c(0, 4),
    check.names = FALSE
  )
  sample_meta <- data.frame(
    Sample = c("S1", "S2", "S3"),
    Group = c("A", "A", NA),
    row.names = c("S1", "S2", "S3"),
    check.names = FALSE
  )

  result <- remove_low_count_genes(
    counts_dat = df,
    sample_metadata = sample_meta,
    feature_id_colname = "Gene",
    group_colname = "Group",
    use_cpm_counts_to_filter = FALSE,
    use_group_based_filtering = TRUE,
    minimum_count_value_to_be_considered_nonzero = 4,
    minimum_number_of_samples_with_nonzero_counts_in_total = 1,
    minimum_number_of_samples_with_nonzero_counts_in_a_group = 2
  )

  expect_equal(result$Gene, "keep")
})

test_that("remove_low_count_genes works with group-based filtering (no grouped tibble crash)", {
  df <- data.frame(
    Gene = c(
      "mt-Nd5_43275",
      "mt-Nd6_43276",
      "mt-Te_43277",
      "mt-Cytb_43278",
      "mt-Tt_43279",
      "mt-Tp_43280"
    ),
    A1 = c(6155, 858, 0, 20542, 0, 12),
    A2 = c(10823, 1420, 1, 29677, 9, 16),
    A3 = c(9482, 1167, 2, 31730, 0, 13),
    B1 = c(6162, 1181, 0, 28293, 0, 15),
    B2 = c(8002, 845, 1, 25617, 7, 19),
    B3 = c(7225, 1198, 3, 30370, 3, 26),
    C1 = c(4141, 515, 0, 21310, 0, 32),
    C2 = c(9058, 1147, 4, 30108, 0, 33),
    C3 = c(8481, 1124, 2, 30893, 2, 50),
    row.names = seq(5, 10)
  )
  sample_meta <- structure(
    list(
      Sample = c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"),
      Group = c("A", "A", "A", "B", "B", "B", "C", "C", "C"),
      Replicate = c(1, 2, 3, 1, 2, 3, 1, 2, 3),
      Batch = c(1, 2, 2, 1, 1, 2, 1, 2, 2),
      Label = c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")
    ),
    row.names = c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"),
    class = "data.frame"
  )

  result <- remove_low_count_genes(
    counts_dat = df,
    sample_metadata = sample_meta,
    feature_id_colname = "Gene",
    group_colname = "Group",
    use_cpm_counts_to_filter = TRUE,
    use_group_based_filtering = TRUE,
    minimum_count_value_to_be_considered_nonzero = 8,
    minimum_number_of_samples_with_nonzero_counts_in_total = 7,
    minimum_number_of_samples_with_nonzero_counts_in_a_group = 3
  )
  expect_s3_class(result, "data.frame")
  expect_true("Gene" %in% colnames(result))
  expect_true(nrow(result) > 0)
})

test_that("filter_counts forwards plotting parameters", {
  pca_args <- NULL
  histogram_args <- NULL

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

  moo <- create_multiOmicDataSet_from_dataframes(
    as.data.frame(nidap_sample_metadata),
    as.data.frame(nidap_clean_raw_counts),
    sample_id_colname = "Sample",
    feature_id_colname = "Gene"
  ) |>
    calc_cpm(feature_id_colname = "Gene")

  filter_counts(
    moo,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    label_colname = "Label",
    count_type = "raw",
    samples_to_rename = c("A1:Alpha 1"),
    add_label_to_pca = FALSE,
    principal_component_on_x_axis = 2,
    principal_component_on_y_axis = 3,
    legend_position_for_pca = "bottom",
    label_offset_x_ = 4,
    label_offset_y_ = 5,
    label_font_size = 6,
    point_size_for_pca = 7,
    color_histogram_by_group = TRUE,
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

  expect_true(histogram_args$color_by_group)
  expect_true(histogram_args$set_min_max_for_x_axis)
  expect_equal(histogram_args$minimum_for_x_axis, -2)
  expect_equal(histogram_args$maximum_for_x_axis, 2)
  expect_equal(histogram_args$x_axis_label, "CPM")
  expect_equal(histogram_args$legend_font_size, 11)
  expect_equal(histogram_args$legend_position, "right")
  expect_equal(histogram_args$number_of_legend_columns, 2)
  expect_equal(
    histogram_args$color_values,
    c(A = "red", B = "blue", C = "green")
  )

  pca_args <- NULL
  filter_counts(
    moo,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    label_colname = "Label",
    count_type = "raw",
    add_label_to_pca = TRUE,
    plot_corr_matrix_heatmap = FALSE,
    print_plots = TRUE,
    save_plots = FALSE
  )
  expect_equal(pca_args$label_colname, "Label")
})

test_that("filter_counts handles histogram label combinations", {
  pca_args <- NULL
  histogram_args <- NULL
  group_colors <- c(A = "red", B = "blue", C = "green")

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

  moo <- create_multiOmicDataSet_from_dataframes(
    as.data.frame(nidap_sample_metadata),
    as.data.frame(nidap_clean_raw_counts),
    sample_id_colname = "Sample",
    feature_id_colname = "Gene"
  ) |>
    calc_cpm(feature_id_colname = "Gene")

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
    filter_counts(
      moo,
      sample_id_colname = "Sample",
      feature_id_colname = "Gene",
      label_colname = combination$label_colname,
      count_type = "raw",
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

test_that("filter_counts forwards the default MOSuite plot colors", {
  pca_args <- NULL
  histogram_args <- NULL
  default_colors <- c(
    "A" = "#5954d6",
    "B" = "#e1562c",
    "C" = "#b80058"
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

  moo <- create_multiOmicDataSet_from_dataframes(
    as.data.frame(nidap_sample_metadata),
    as.data.frame(nidap_clean_raw_counts),
    sample_id_colname = "Sample",
    feature_id_colname = "Gene"
  ) |>
    calc_cpm(feature_id_colname = "Gene")

  filter_counts(
    moo,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    label_colname = "Label",
    count_type = "raw",
    plot_corr_matrix_heatmap = FALSE,
    print_plots = TRUE,
    save_plots = FALSE
  )

  expect_equal(pca_args$color_values, default_colors)
  expect_equal(histogram_args$color_values, default_colors)
})

test_that("filter_counts adds a dotted threshold line and label to the histogram", {
  histogram_plot <- NULL
  histogram_args <- NULL

  local_mocked_bindings(
    plot_pca = function(...) {
      return(ggplot2::ggplot())
    },
    plot_histogram = function(...) {
      histogram_args <<- list(...)
      return(ggplot2::ggplot())
    },
    print_or_save_plot = function(plot, filename, ...) {
      if (basename(filename) == "histogram.png") {
        histogram_plot <<- plot
      }
      return(invisible(NULL))
    },
    .package = "MOSuite"
  )

  moo <- create_multiOmicDataSet_from_dataframes(
    as.data.frame(nidap_sample_metadata),
    as.data.frame(nidap_clean_raw_counts),
    sample_id_colname = "Sample",
    feature_id_colname = "Gene"
  )

  filter_counts(
    moo,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    label_colname = "Label",
    count_type = "raw",
    minimum_count_value_to_be_considered_nonzero = 4,
    plot_corr_matrix_heatmap = FALSE,
    print_plots = TRUE,
    save_plots = FALSE
  )

  histogram_data <- ggplot2::ggplot_build(histogram_plot)$data

  expect_true(histogram_args$use_log2_x_axis)
  expect_true(histogram_args$return_ggplot)
  expect_equal(histogram_data[[1]]$xintercept, 4 + 0.5)
  expect_equal(histogram_data[[1]]$linetype, 2)
  expect_equal(histogram_data[[1]]$linewidth, 1)
  expect_equal(histogram_data[[2]]$x, 4 + 0.5)
  expect_equal(histogram_data[[2]]$label, "CPM: 4")
})

test_that("filter_counts PCA matches standalone plot_pca on filtered output", {
  pca_capture <- capture_saved_pca_plot()
  local_mocked_bindings(
    print_or_save_plot = pca_capture$print_or_save_plot,
    .package = "MOSuite"
  )

  moo <- create_multiOmicDataSet_from_dataframes(
    as.data.frame(nidap_sample_metadata),
    as.data.frame(nidap_clean_raw_counts),
    sample_id_colname = "Sample",
    feature_id_colname = "Gene"
  ) |>
    filter_counts(
      sample_id_colname = "Sample",
      feature_id_colname = "Gene",
      label_colname = NULL,
      count_type = "raw",
      use_cpm_counts_to_filter = FALSE,
      minimum_count_value_to_be_considered_nonzero = 8,
      minimum_number_of_samples_with_nonzero_counts_in_total = 7,
      plot_corr_matrix_heatmap = FALSE,
      print_plots = TRUE,
      save_plots = FALSE
    )

  expected_pca <- plot_pca(
    moo@counts$filt,
    sample_metadata = moo@sample_meta,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    group_colname = "Group",
    label_colname = NULL,
    samples_to_rename = c(""),
    principal_components = c(1, 2),
    legend_position = "top",
    point_size = 5,
    label_font_size = 3,
    label_offset_y_ = 2,
    label_offset_x_ = 2,
    log_transform = TRUE,
    log_transform_pseudocount = 0.5,
    log_transform_base = "ln",
    print_plots = FALSE,
    save_plots = FALSE
  )

  expect_s3_class(pca_capture$get(), "ggplot")
  expect_pca_coordinates_equal(pca_capture$get(), expected_pca)
})

test_that("filter_counts histogram matches standalone plot_histogram on filtered output", {
  histogram_capture <- capture_saved_histogram_plot()
  local_mocked_bindings(
    print_or_save_plot = histogram_capture$print_or_save_plot,
    .package = "MOSuite"
  )

  moo <- create_multiOmicDataSet_from_dataframes(
    as.data.frame(nidap_sample_metadata),
    as.data.frame(nidap_clean_raw_counts),
    sample_id_colname = "Sample",
    feature_id_colname = "Gene"
  ) |>
    filter_counts(
      sample_id_colname = "Sample",
      feature_id_colname = "Gene",
      label_colname = NULL,
      count_type = "raw",
      use_cpm_counts_to_filter = FALSE,
      minimum_count_value_to_be_considered_nonzero = 8,
      minimum_number_of_samples_with_nonzero_counts_in_total = 7,
      plot_corr_matrix_heatmap = FALSE,
      print_plots = TRUE,
      save_plots = FALSE,
      interactive_plots = FALSE
    )

  histogram_threshold_x <- 8 + 0.5
  expected_histogram <- plot_histogram(
    moo@counts$filt,
    sample_metadata = moo@sample_meta,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    group_colname = "Group",
    label_colname = NULL,
    color_by_group = TRUE,
    set_min_max_for_x_axis = FALSE,
    minimum_for_x_axis = -1,
    maximum_for_x_axis = 1,
    x_axis_label = "Count",
    legend_position = "top",
    legend_font_size = NULL,
    number_of_legend_columns = 6,
    interactive_plots = FALSE,
    return_ggplot = TRUE,
    use_log2_x_axis = TRUE
  ) +
    ggplot2::labs(caption = "filtered counts") +
    ggplot2::geom_vline(
      xintercept = histogram_threshold_x,
      linetype = 2,
      linewidth = 1
    ) +
    ggplot2::annotate(
      "text",
      x = histogram_threshold_x,
      y = Inf,
      label = "Count: 8",
      hjust = -0.05,
      vjust = 1.5,
      size = 3
    )

  expect_s3_class(histogram_capture$get(), "ggplot")
  expect_histogram_layers_equal(histogram_capture$get(), expected_histogram)
})
