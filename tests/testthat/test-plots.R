set.seed(20250225)

corr_heatmap_fixture <- function() {
  plot_corr_heatmap(
    nidap_filtered_counts |>
      as.data.frame(),
    sample_metadata = as.data.frame(nidap_sample_metadata),
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    label_colname = "Label",
    group_colname = "Group",
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
    )
  )
}

test_that("print_or_save_plot saves ComplexHeatmap to disk without error", {
  p <- corr_heatmap_fixture()
  outfile <- tempfile(fileext = ".png")
  result <- print_or_save_plot(
    p,
    filename = outfile,
    print_plots = FALSE,
    save_plots = TRUE,
    plots_dir = "",
    caption = "filtered counts"
  )
  expect_equal(result, outfile)
  expect_true(file.exists(outfile))
  expect_gt(file.size(outfile), 0)
})

test_that("print_or_save_plot saves ggplot without error", {
  p <- plot_read_depth(nidap_clean_raw_counts)
  outfile <- tempfile(fileext = ".png")
  result <- print_or_save_plot(
    p,
    filename = outfile,
    print_plots = FALSE,
    save_plots = TRUE,
    plots_dir = "",
    caption = "normalized counts"
  )
  expect_equal(result, outfile)
  expect_true(file.exists(outfile))
  expect_gt(file.size(outfile), 0)
})

test_that("print_or_save_plot prints ComplexHeatmap with caption without error", {
  p <- corr_heatmap_fixture()
  outfile <- tempfile(fileext = ".png")
  withr::with_png(outfile, {
    result <- print_or_save_plot(
      p,
      filename = outfile,
      print_plots = TRUE,
      save_plots = FALSE,
      plots_dir = "",
      caption = "batch-corrected counts"
    )
  })
  expect_equal(result, outfile)
})

test_that("save_or_print_plot works for ComplexHeatmap", {
  p <- corr_heatmap_fixture()
  skip_on_ci()
  expect_snapshot_file(
    print_or_save_plot(
      p,
      filename = "heatmap.png",
      print_plots = FALSE,
      save_plots = TRUE,
      plots_dir = "."
    ),
    "heatmap.png"
  )
})
test_that("save_or_print_plot works for ggplot", {
  p <- plot_read_depth(nidap_clean_raw_counts)
  skip_on_ci()
  expect_snapshot_file(
    print_or_save_plot(
      p,
      filename = "read_depth.png",
      print_plots = FALSE,
      save_plots = TRUE,
      plots_dir = "."
    ),
    "read_depth.png"
  )
})
