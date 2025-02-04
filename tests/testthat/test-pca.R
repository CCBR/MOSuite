test_that("calc_pca works", {
  pca_dat <- calc_pca(
    nidap_clean_raw_counts,
    nidap_sample_metadata
  ) %>%
    dplyr::filter(PC %in% c(1, 2))
  expect_equal(
    pca_dat,
    structure(list(Sample = c(
      "A1", "A1", "A2", "A2", "A3", "A3",
      "B1", "B1", "B2", "B2", "B3", "B3", "C1", "C1", "C2", "C2", "C3",
      "C3"
    ), PC = c(
      1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
      1, 2
    ), value = c(
      -40.6241668816455, 25.2297268619146, -56.2133160433603,
      6.13385771612248, -69.1070711020441, -21.8952345106934, -36.1660251215743,
      7.80504297978752, -25.865255255388, -11.2138080494717, -9.6232450176941,
      9.32724696042314, 74.3345576680281, -86.7286802229905, 85.0442226808989,
      117.992340543509, 78.2202990727852, -46.6504922786012
    ), std.dev = c(
      61.7780383925471,
      55.9548424792563, 61.7780383925471, 55.9548424792563, 61.7780383925471,
      55.9548424792563, 61.7780383925471, 55.9548424792563, 61.7780383925471,
      55.9548424792563, 61.7780383925471, 55.9548424792563, 61.7780383925471,
      55.9548424792563, 61.7780383925471, 55.9548424792563, 61.7780383925471,
      55.9548424792563
    ), percent = c(
      21.219, 17.408, 21.219, 17.408,
      21.219, 17.408, 21.219, 17.408, 21.219, 17.408, 21.219, 17.408,
      21.219, 17.408, 21.219, 17.408, 21.219, 17.408
    ), cumulative = c(
      0.21219,
      0.38627, 0.21219, 0.38627, 0.21219, 0.38627, 0.21219, 0.38627,
      0.21219, 0.38627, 0.21219, 0.38627, 0.21219, 0.38627, 0.21219,
      0.38627, 0.21219, 0.38627
    ), Group = c(
      "A", "A", "A", "A", "A",
      "A", "B", "B", "B", "B", "B", "B", "C", "C", "C", "C", "C", "C"
    ), Replicate = c(
      1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3, 1, 1, 2,
      2, 3, 3
    ), Batch = c(
      1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1,
      2, 2, 2, 2
    ), Label = c(
      "A1", "A1", "A2", "A2", "A3", "A3", "B1",
      "B1", "B2", "B2", "B3", "B3", "C1", "C1", "C2", "C2", "C3", "C3"
    )), class = c("tbl_df", "tbl", "data.frame"), row.names = c(
      NA,
      -18L
    ))
  )
})

test_that("plot_pca layers are expected", {
  p <- plot_pca(
    counts_dat = nidap_filtered_counts,
    sample_metadata = nidap_sample_metadata,
    samples_to_rename = NULL,
    group_colname = "Group",
    label_colname = "Label",
    color_values = c(
      "#5954d6", "#e1562c", "#b80058", "#00c6f8", "#d163e6", "#00a76c",
      "#ff9287", "#008cf9", "#006e00", "#796880", "#FFA500", "#878500"
    ),
    principal_components = c(1, 2),
    legend_position_for_pca = "top",
    point_size_for_pca = 1,
    add_label_to_pca = TRUE,
    label_font_size = 3,
    label_offset_y_ = 2,
    label_offset_x_ = 2
  )

  expect_s3_class(p$layers[[1]], "ggproto")
  expect_s3_class(p$layers[[1]]$geom, "GeomPoint")
})

test_that("3DPCA works", {
  skip()
  plot_pca_3d(nidap_filtered_counts,
    nidap_sample_metadata,
    group_colname = "Group",
    label_colname = "Label",
    color_values = c(
      "#5954d6", "#e1562c", "#b80058", "#00c6f8", "#d163e6", "#00a76c",
      "#ff9287", "#008cf9", "#006e00", "#796880", "#FFA500", "#878500"
    ),
    principal_components = c(1, 2, 3),
  )
})
