test_that("calc_pca works", {
  pca_dat <- calc_pca(nidap_clean_raw_counts, nidap_sample_metadata) |>
    dplyr::filter(PC %in% c(1, 2))
  expect_equal(
    pca_dat,
    structure(
      list(
        Sample = c(
          "A1",
          "A1",
          "A2",
          "A2",
          "A3",
          "A3",
          "B1",
          "B1",
          "B2",
          "B2",
          "B3",
          "B3",
          "C1",
          "C1",
          "C2",
          "C2",
          "C3",
          "C3"
        ),
        PC = c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2),
        value = c(
          -40.6241668816455,
          25.2297268619146,
          -56.2133160433603,
          6.13385771612248,
          -69.1070711020441,
          -21.8952345106934,
          -36.1660251215743,
          7.80504297978752,
          -25.865255255388,
          -11.2138080494717,
          -9.6232450176941,
          9.32724696042314,
          74.3345576680281,
          -86.7286802229905,
          85.0442226808989,
          117.992340543509,
          78.2202990727852,
          -46.6504922786012
        ),
        std.dev = c(
          61.7780383925471,
          55.9548424792563,
          61.7780383925471,
          55.9548424792563,
          61.7780383925471,
          55.9548424792563,
          61.7780383925471,
          55.9548424792563,
          61.7780383925471,
          55.9548424792563,
          61.7780383925471,
          55.9548424792563,
          61.7780383925471,
          55.9548424792563,
          61.7780383925471,
          55.9548424792563,
          61.7780383925471,
          55.9548424792563
        ),
        percent = c(
          21.219,
          17.408,
          21.219,
          17.408,
          21.219,
          17.408,
          21.219,
          17.408,
          21.219,
          17.408,
          21.219,
          17.408,
          21.219,
          17.408,
          21.219,
          17.408,
          21.219,
          17.408
        ),
        cumulative = c(
          0.21219,
          0.38627,
          0.21219,
          0.38627,
          0.21219,
          0.38627,
          0.21219,
          0.38627,
          0.21219,
          0.38627,
          0.21219,
          0.38627,
          0.21219,
          0.38627,
          0.21219,
          0.38627,
          0.21219,
          0.38627
        ),
        Group = c(
          "A",
          "A",
          "A",
          "A",
          "A",
          "A",
          "B",
          "B",
          "B",
          "B",
          "B",
          "B",
          "C",
          "C",
          "C",
          "C",
          "C",
          "C"
        ),
        Replicate = c(1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3),
        Batch = c(1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 2, 2),
        Label = c(
          "A1",
          "A1",
          "A2",
          "A2",
          "A3",
          "A3",
          "B1",
          "B1",
          "B2",
          "B2",
          "B3",
          "B3",
          "C1",
          "C1",
          "C2",
          "C2",
          "C3",
          "C3"
        )
      ),
      class = c("tbl_df", "tbl", "data.frame"),
      row.names = c(NA, -18L)
    )
  )
})

test_that("plot_pca layers are expected", {
  p <- plot_pca(
    moo_counts = nidap_filtered_counts,
    sample_metadata = nidap_sample_metadata,
    principal_components = c(1, 2),
    samples_to_rename = NULL,
    group_colname = "Group",
    label_colname = "Label",
    color_values = c(
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
    ),
    legend_position = "top",
    point_size = 1,
    add_label = TRUE,
    label_font_size = 3,
    label_offset_y_ = 2,
    label_offset_x_ = 2
  )

  expect_s3_class(p$layers[[1]], "ggproto")
  expect_s3_class(p$layers[[1]]$geom, "GeomPoint")
})

normalize_color_values <- function(colors) {
  vapply(
    colors,
    function(color) {
      if (grepl("^rgba\\(", color)) {
        color_parts <- strsplit(gsub("^rgba\\(|\\)$", "", color), ",")[[1]]
        color_parts <- as.numeric(color_parts[seq_len(3)])
        grDevices::rgb(
          color_parts[1],
          color_parts[2],
          color_parts[3],
          maxColorValue = 255
        )
      } else {
        rgb_value <- grDevices::col2rgb(color)
        grDevices::rgb(
          rgb_value[1, 1],
          rgb_value[2, 1],
          rgb_value[3, 1],
          maxColorValue = 255
        )
      }
    },
    character(1),
    USE.NAMES = FALSE
  )
}

get_colour_scale <- function(plot) {
  scales <- ggplot2::ggplot_build(plot)$plot$scales$scales
  scales[[which(vapply(
    scales,
    function(scale) "colour" %in% scale$aesthetics,
    logical(1)
  ))[[1]]]]
}

get_colour_guide_ncol <- function(plot) {
  plot$guides$guides$colour$params$ncol
}

test_that("2D PCA wraps long top and bottom sample-name legends", {
  sample_columns <- setdiff(colnames(nidap_filtered_counts), "Gene")
  long_sample_names <- stats::setNames(
    sprintf("SampleName%05d", seq_along(sample_columns)),
    sample_columns
  )
  counts_dat <- nidap_filtered_counts
  colnames(counts_dat) <- ifelse(
    colnames(counts_dat) %in% names(long_sample_names),
    unname(long_sample_names[colnames(counts_dat)]),
    colnames(counts_dat)
  )
  sample_metadata <- nidap_sample_metadata
  sample_metadata$Sample <- unname(long_sample_names[as.character(
    sample_metadata$Sample
  )])
  sample_metadata$Label <- sample_metadata$Sample

  for (legend_position in c("top", "bottom")) {
    pca_2d <- plot_pca_2d(
      counts_dat,
      sample_metadata = sample_metadata,
      sample_id_colname = "Sample",
      feature_id_colname = "Gene",
      group_colname = "Sample",
      label_colname = "Label",
      legend_position = legend_position,
      add_label = FALSE,
      print_plots = FALSE,
      save_plots = FALSE
    )

    expect_equal(get_colour_guide_ncol(pca_2d), 3)
  }
})

test_that("2D and 3D PCA resolve unnamed colors by first observed group order", {
  color_values <- c("#5954d6", "#e1562c", "#b80058")
  expected_colors <- c(B = "#5954d6", A = "#e1562c", C = "#b80058")
  counts_dat <- nidap_filtered_counts[, c(
    "Gene",
    "B1",
    "B2",
    "B3",
    "A1",
    "A2",
    "A3",
    "C1",
    "C2",
    "C3"
  )]

  pca_2d <- plot_pca_2d(
    counts_dat,
    sample_metadata = nidap_sample_metadata,
    feature_id_colname = "Gene",
    color_values = color_values,
    add_label = FALSE,
    print_plots = FALSE,
    save_plots = FALSE
  )
  pca_2d_colors <- get_colour_scale(pca_2d)$palette.cache[names(
    expected_colors
  )]

  pca_3d <- plot_pca_3d(
    counts_dat,
    sample_metadata = nidap_sample_metadata,
    feature_id_colname = "Gene",
    color_values = color_values,
    print_plots = FALSE,
    save_plots = FALSE
  )
  pca_3d_traces <- plotly::plotly_build(pca_3d)$x$data
  pca_3d_colors <- stats::setNames(
    normalize_color_values(vapply(
      pca_3d_traces,
      function(trace) trace$marker$color,
      character(1)
    )),
    vapply(pca_3d_traces, function(trace) trace$name, character(1))
  )[names(expected_colors)]

  expected_colors <- stats::setNames(
    normalize_color_values(expected_colors),
    names(expected_colors)
  )

  expect_equal(normalize_color_values(pca_2d_colors), unname(expected_colors))
  expect_equal(pca_3d_colors, expected_colors)
})

test_that("2D and 3D PCA resolve unnamed colors by factor level order", {
  color_values <- c("#5954d6", "#e1562c", "#b80058")
  expected_colors <- c(C = "#5954d6", A = "#e1562c", B = "#b80058")
  counts_dat <- nidap_filtered_counts[, c(
    "Gene",
    "B1",
    "B2",
    "B3",
    "A1",
    "A2",
    "A3",
    "C1",
    "C2",
    "C3"
  )]
  sample_metadata <- nidap_sample_metadata
  sample_metadata$Group <- factor(
    sample_metadata$Group,
    levels = c("C", "A", "B")
  )

  pca_2d <- plot_pca_2d(
    counts_dat,
    sample_metadata = sample_metadata,
    feature_id_colname = "Gene",
    color_values = color_values,
    add_label = FALSE,
    print_plots = FALSE,
    save_plots = FALSE
  )
  pca_2d_colors <- get_colour_scale(pca_2d)$palette.cache[names(
    expected_colors
  )]

  pca_3d <- plot_pca_3d(
    counts_dat,
    sample_metadata = sample_metadata,
    feature_id_colname = "Gene",
    color_values = color_values,
    print_plots = FALSE,
    save_plots = FALSE
  )
  pca_3d_traces <- plotly::plotly_build(pca_3d)$x$data
  pca_3d_colors <- stats::setNames(
    normalize_color_values(vapply(
      pca_3d_traces,
      function(trace) trace$marker$color,
      character(1)
    )),
    vapply(pca_3d_traces, function(trace) trace$name, character(1))
  )[names(expected_colors)]

  expected_colors <- stats::setNames(
    normalize_color_values(expected_colors),
    names(expected_colors)
  )

  expect_equal(normalize_color_values(pca_2d_colors), unname(expected_colors))
  expect_equal(pca_3d_colors, expected_colors)
})

test_that("2D PCA preserves named color mappings", {
  color_values <- c(C = "#5954d6", A = "#e1562c", B = "#b80058")

  pca_2d <- plot_pca_2d(
    nidap_filtered_counts,
    sample_metadata = nidap_sample_metadata,
    feature_id_colname = "Gene",
    color_values = color_values,
    add_label = FALSE,
    print_plots = FALSE,
    save_plots = FALSE
  )

  expect_equal(get_colour_scale(pca_2d)$palette.cache, color_values)
})


test_that("2D & 3D PCA method dispatch works", {
  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts)
    )
  )
  expect_equal(
    plot_pca(
      moo,
      count_type = "filt",
      principal_components = c(1, 2)
    ),
    plot_pca(
      moo@counts$filt,
      moo@sample_meta,
      principal_components = c(1, 2)
    )
  )

  # 3D PCA
  p1 <- plot_pca(moo, count_type = "filt", principal_components = c(1, 2, 3))
  p2 <- plot_pca(
    moo@counts$filt,
    moo@sample_meta,
    principal_components = c(1, 2, 3)
  )
  # see compare_proxy.plotly
  # expect_equal(p1, p2)
})

test_that("plot_pca_3d returns plotly object and has correct structure", {
  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts)
    )
  )

  # Test with multiOmicDataSet
  fig_moo <- plot_pca_3d(
    moo,
    count_type = "filt",
    principal_components = c(1, 2, 3),
    group_colname = "Group",
    label_colname = "Label",
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(fig_moo, "plotly")
  expect_type(fig_moo$x, "list")

  # Test with data.frame
  fig_df <- plot_pca_3d(
    moo@counts$filt,
    sample_metadata = moo@sample_meta,
    principal_components = c(1, 2, 3),
    group_colname = "Group",
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(fig_df, "plotly")
  expect_type(fig_df$x, "list")
})

test_that("plot_pca_3d validates principal_components length", {
  expect_error(
    plot_pca_3d(
      nidap_filtered_counts,
      sample_metadata = nidap_sample_metadata,
      principal_components = c(1, 2),
      save_plots = FALSE,
      print_plots = FALSE
    ),
    "principal_components must contain 3 values"
  )

  expect_error(
    plot_pca_3d(
      nidap_filtered_counts,
      sample_metadata = nidap_sample_metadata,
      principal_components = c(1, 2, 3, 4),
      save_plots = FALSE,
      print_plots = FALSE
    ),
    "principal_components must contain 3 values"
  )
})

test_that("plot_pca_2d works on multiOmicDataSet object", {
  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts)
    )
  )

  # Test with multiOmicDataSet
  p_moo <- plot_pca_2d(
    moo,
    count_type = "filt",
    principal_components = c(1, 2),
    group_colname = "Group",
    label_colname = "Label",
    save_plots = FALSE,
    print_plots = FALSE
  )

  expect_s3_class(p_moo, "ggplot")
  # Should have geom_point and geom_text_repel layers
  expect_gte(length(p_moo$layers), 2)
  expect_s3_class(p_moo$layers[[1]]$geom, "GeomPoint")
})

test_that("plot_pca_2d works with and without labels", {
  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts)
    )
  )

  # With labels
  p_with_labels <- plot_pca_2d(
    moo,
    count_type = "filt",
    principal_components = c(1, 2),
    add_label = TRUE,
    save_plots = FALSE,
    print_plots = FALSE
  )

  # Without labels
  p_without_labels <- plot_pca_2d(
    moo,
    count_type = "filt",
    principal_components = c(1, 2),
    add_label = FALSE,
    save_plots = FALSE,
    print_plots = FALSE
  )

  # With labels should have more layers (geom_text_repel)
  expect_gt(length(p_with_labels$layers), length(p_without_labels$layers))
})
