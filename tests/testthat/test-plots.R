test_that("save_or_print_plot works for ComplexHeatmap", {
  p <- plot_corr_heatmap(
    counts_dat = nidap_filtered_counts %>%
      dplyr::select(tidyselect::all_of(
        c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")
      )) %>%
      as.data.frame(),
    sample_metadata = as.data.frame(nidap_sample_metadata),
    sample_id_colname = "Sample",
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
  expect_snapshot_file(
    print_or_save_plot(p,
      filename = "heatmap.png",
      print_plots = FALSE, save_plots = TRUE, plots_dir = ""
    ),
    "heatmap.png"
  )
})
test_that("save_or_print_plot works for ggplot", {
  p <- plot_read_depth(nidap_clean_raw_counts)
  expect_snapshot_file(
    print_or_save_plot(p,
      filename = "read_depth.png",
      print_plots = TRUE, save_plots = TRUE, plots_dir = NULL
    ),
    "read_depth.png"
  )
})
