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
