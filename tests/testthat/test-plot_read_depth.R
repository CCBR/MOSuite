test_that("plot_read_depth works on moo & dataframes", {
  moo <- multiOmicDataSet(
    sample_metadata = nidap_sample_metadata,
    anno_dat = data.frame(),
    counts_lst = list("raw" = nidap_raw_counts)
  )
  expect_equal(
    plot_read_depth(moo, "raw"),
    plot_read_depth(nidap_raw_counts)
  )
})

test_that("plot_read_depth accepts extra args via moo dispatch without error", {
  moo <- multiOmicDataSet(
    sample_metadata = nidap_sample_metadata,
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = nidap_raw_counts,
      "clean" = nidap_clean_raw_counts
    )
  )
  expect_no_error(
    plot_read_depth(
      moo,
      count_type = "clean",
      sample_id_colname = "Sample",
      feature_id_colname = "Gene"
    )
  )
})

test_that("plot_read_depth can color samples by group", {
  moo <- multiOmicDataSet(
    sample_metadata = nidap_sample_metadata,
    anno_dat = data.frame(),
    counts_lst = list("raw" = nidap_raw_counts)
  )
  plot <- plot_read_depth(
    moo,
    count_type = "raw",
    sample_id_colname = "Sample",
    group_colname = "Group",
    color_values = c("A" = "blue", "B" = "green", "C" = "orange")
  )

  built <- ggplot2::ggplot_build(plot)
  expect_equal(
    unique(built$data[[1]]$fill),
    c("blue", "green", "orange")
  )
})

test_that("plot_read_depth keeps single-color bars when group_colname is blank", {
  moo <- multiOmicDataSet(
    sample_metadata = nidap_sample_metadata,
    anno_dat = data.frame(),
    counts_lst = list("raw" = nidap_raw_counts)
  )
  plot <- plot_read_depth(
    moo,
    count_type = "raw",
    sample_id_colname = "Sample",
    group_colname = "",
    color_values = c("A" = "blue", "B" = "green", "C" = "orange")
  )

  built <- ggplot2::ggplot_build(plot)
  expect_equal(unique(built$data[[1]]$fill), "blue")
})

test_that("plot_read_depth extends undersupplied group colors", {
  moo <- multiOmicDataSet(
    sample_metadata = nidap_sample_metadata,
    anno_dat = data.frame(),
    counts_lst = list("raw" = nidap_raw_counts)
  )
  plot <- plot_read_depth(
    moo,
    count_type = "raw",
    sample_id_colname = "Sample",
    group_colname = "Group",
    color_values = c("blue")
  )

  built <- ggplot2::ggplot_build(plot)
  expect_equal(
    unique(built$data[[1]]$fill),
    c("blue", "#e1562c", "#b80058")
  )
})
