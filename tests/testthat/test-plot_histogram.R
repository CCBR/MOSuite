log_counts <- structure(
  c(
    4.68203758078952,
    4.85622577980595,
    4.78525385564269,
    4.49631900610197,
    4.06045359796882,
    4.23984554358195,
    4.5945977606547,
    4.64452454872198,
    5.01435137189661,
    4.97202092328668,
    4.63530665635318,
    4.2685504137767,
    4.55963475259408,
    4.48551456875509,
    4.8871583679813,
    4.80572893488115,
    4.62485692334614,
    3.28500783834613,
    4.27910086185478,
    4.44991091582054,
    5.1543875974913,
    5.04097648388479,
    3.93341992672261,
    3.60201340344135,
    4.38819757498796,
    4.52837375303351,
    5.04374399504142,
    5.16774070067169,
    4.88762227733767,
    3.4002608223214,
    4.86037842060906,
    5.24497940734463,
    4.91946017962122,
    5.04632856830695,
    4.96825190546972,
    4.23362208987531,
    4.11790170891073,
    4.05303076635639,
    4.92868500786216,
    4.94726383191083,
    3.80832003906478,
    4.59490650795107,
    4.85738793032812,
    4.40655427188147,
    4.22382743940071,
    4.58770826515446,
    4.96057466313368,
    4.99824283163569,
    4.5437556841171,
    4.97124300205374,
    4.93192370480255,
    4.13835014251084,
    5.0939015418123,
    5.03285949305828
  ),
  dim = c(6L, 9L),
  dimnames = list(
    c(
      "Mrpl15_32",
      "Lypla1_34",
      "Tcea1_36",
      "Atp6v1h_44",
      "Rb1cc1_54",
      "Pcmtd1_68"
    ),
    c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")
  )
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

get_histogram_colour_guide_ncol <- function(plot) {
  return(plot$guides$guides$colour$params$ncol)
}

get_plotly_text <- function(plot) {
  traces <- plotly::plotly_build(plot)$x$data
  return(unlist(
    lapply(traces, function(trace) trace$text),
    use.names = FALSE
  ))
}

test_that("plot_histogram interactive hover text includes sample and group", {
  counts_dat <- log_counts |>
    as.data.frame() |>
    tibble::rownames_to_column("Gene")

  plot_by_group <- suppressWarnings(plot_histogram(
    counts_dat,
    sample_metadata = sample_meta,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    group_colname = "Group",
    color_by_group = TRUE,
    interactive_plots = TRUE
  ))
  plot_by_sample <- suppressWarnings(plot_histogram(
    counts_dat,
    sample_metadata = sample_meta,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    group_colname = "Group",
    color_by_group = FALSE,
    interactive_plots = TRUE
  ))

  hover_text_by_group <- get_plotly_text(plot_by_group)
  hover_text_by_sample <- get_plotly_text(plot_by_sample)
  expect_true(any(grepl("Sample: A1", hover_text_by_group, fixed = TRUE)))
  expect_true(any(grepl("Group: A", hover_text_by_group, fixed = TRUE)))
  expect_true(any(grepl("Sample: A1", hover_text_by_sample, fixed = TRUE)))
  expect_true(any(grepl("Group: A", hover_text_by_sample, fixed = TRUE)))
})

test_that("plot_histogram interactive output keeps legend", {
  counts_dat <- log_counts |>
    as.data.frame() |>
    tibble::rownames_to_column("Gene")

  plot <- suppressWarnings(plot_histogram(
    counts_dat,
    sample_metadata = sample_meta,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    group_colname = "Group",
    color_by_group = TRUE,
    interactive_plots = TRUE
  ))
  built_plot <- plotly::plotly_build(plot)

  expect_false(identical(built_plot$x$layout$showlegend, FALSE))
  expect_true(any(vapply(
    built_plot$x$data,
    function(trace) !is.null(trace$name) && nzchar(trace$name),
    logical(1)
  )))
})

test_that("plot_histogram wraps long top and bottom sample-name legends", {
  counts_dat <- log_counts |>
    as.data.frame() |>
    tibble::rownames_to_column("Gene")
  sample_columns <- setdiff(colnames(counts_dat), "Gene")
  long_sample_names <- stats::setNames(
    sprintf("SampleName%05d", seq_along(sample_columns)),
    sample_columns
  )
  colnames(counts_dat) <- ifelse(
    colnames(counts_dat) %in% names(long_sample_names),
    unname(long_sample_names[colnames(counts_dat)]),
    colnames(counts_dat)
  )
  sample_metadata <- sample_meta
  sample_metadata$Sample <- unname(long_sample_names[as.character(
    sample_metadata$Sample
  )])
  sample_metadata$Label <- sample_metadata$Sample

  for (legend_position in c("top", "bottom")) {
    plot <- plot_histogram(
      counts_dat,
      sample_metadata = sample_metadata,
      sample_id_colname = "Sample",
      feature_id_colname = "Gene",
      color_by_group = FALSE,
      legend_position = legend_position,
      number_of_legend_columns = 6
    )

    expect_equal(get_histogram_colour_guide_ncol(plot), 3)
  }
})

test_that("plot_histogram legend columns target the colour guide", {
  plot <- plot_histogram(
    log_counts |>
      as.data.frame() |>
      tibble::rownames_to_column("Gene"),
    sample_meta,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    color_by_group = FALSE,
    legend_position = "top",
    number_of_legend_columns = 2
  )

  expect_equal(get_histogram_colour_guide_ncol(plot), 2)
})

test_that("plot_histogram uses line glyphs for density legend keys", {
  plot <- plot_histogram(
    nidap_filtered_counts,
    sample_metadata = nidap_sample_metadata,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    group_colname = "Group",
    color_by_group = TRUE
  )

  expect_equal(
    plot$guides$guides$colour$params$override.aes,
    list(linetype = 1, linewidth = 2, shape = NA, fill = NA)
  )
  expect_s3_class(plot$theme$legend.key, "element_blank")
})

test_that("plot_histogram works with rownames", {
  p <- plot_histogram(
    log_counts |> as.data.frame() |> tibble::rownames_to_column("Gene"),
    sample_meta,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    group_colname = "Group",
    label_colname = "Label",
    color_values = c(
      indigo = "#5954d6",
      carrot = "#e1562c",
      lipstick = "#b80058",
      turquoise = "#00c6f8",
      lavender = "#d163e6",
      jade = "#00a76c",
      coral = "#ff9287",
      azure = "#008cf9",
      green = "#006e00",
      rum = "#796880",
      orange = "#FFA500",
      olive = "#878500"
    ),
    color_by_group = FALSE,
    set_min_max_for_x_axis = FALSE,
    minimum_for_x_axis = -1,
    maximum_for_x_axis = 1,
    legend_position = "top",
    legend_font_size = 10,
    number_of_legend_columns = 6
  )

  expect_s3_class(p$layers[[1]], "ggproto")
  expect_s3_class(p$layers[[1]]$geom, "GeomArea")
  expect_equal(
    p$data,
    structure(
      list(
        Gene = c(
          "Mrpl15_32",
          "Mrpl15_32",
          "Mrpl15_32",
          "Mrpl15_32",
          "Mrpl15_32",
          "Mrpl15_32",
          "Mrpl15_32",
          "Mrpl15_32",
          "Mrpl15_32",
          "Lypla1_34",
          "Lypla1_34",
          "Lypla1_34",
          "Lypla1_34",
          "Lypla1_34",
          "Lypla1_34",
          "Lypla1_34",
          "Lypla1_34",
          "Lypla1_34",
          "Tcea1_36",
          "Tcea1_36",
          "Tcea1_36",
          "Tcea1_36",
          "Tcea1_36",
          "Tcea1_36",
          "Tcea1_36",
          "Tcea1_36",
          "Tcea1_36",
          "Atp6v1h_44",
          "Atp6v1h_44",
          "Atp6v1h_44",
          "Atp6v1h_44",
          "Atp6v1h_44",
          "Atp6v1h_44",
          "Atp6v1h_44",
          "Atp6v1h_44",
          "Atp6v1h_44",
          "Rb1cc1_54",
          "Rb1cc1_54",
          "Rb1cc1_54",
          "Rb1cc1_54",
          "Rb1cc1_54",
          "Rb1cc1_54",
          "Rb1cc1_54",
          "Rb1cc1_54",
          "Rb1cc1_54",
          "Pcmtd1_68",
          "Pcmtd1_68",
          "Pcmtd1_68",
          "Pcmtd1_68",
          "Pcmtd1_68",
          "Pcmtd1_68",
          "Pcmtd1_68",
          "Pcmtd1_68",
          "Pcmtd1_68"
        ),
        Sample = c(
          "A1",
          "A2",
          "A3",
          "B1",
          "B2",
          "B3",
          "C1",
          "C2",
          "C3",
          "A1",
          "A2",
          "A3",
          "B1",
          "B2",
          "B3",
          "C1",
          "C2",
          "C3",
          "A1",
          "A2",
          "A3",
          "B1",
          "B2",
          "B3",
          "C1",
          "C2",
          "C3",
          "A1",
          "A2",
          "A3",
          "B1",
          "B2",
          "B3",
          "C1",
          "C2",
          "C3",
          "A1",
          "A2",
          "A3",
          "B1",
          "B2",
          "B3",
          "C1",
          "C2",
          "C3",
          "A1",
          "A2",
          "A3",
          "B1",
          "B2",
          "B3",
          "C1",
          "C2",
          "C3"
        ),
        count = c(
          4.68203758078952,
          4.5945977606547,
          4.55963475259408,
          4.27910086185478,
          4.38819757498796,
          4.86037842060906,
          4.11790170891073,
          4.85738793032812,
          4.5437556841171,
          4.85622577980595,
          4.64452454872198,
          4.48551456875509,
          4.44991091582054,
          4.52837375303351,
          5.24497940734463,
          4.05303076635639,
          4.40655427188147,
          4.97124300205374,
          4.78525385564269,
          5.01435137189661,
          4.8871583679813,
          5.1543875974913,
          5.04374399504142,
          4.91946017962122,
          4.92868500786216,
          4.22382743940071,
          4.93192370480255,
          4.49631900610197,
          4.97202092328668,
          4.80572893488115,
          5.04097648388479,
          5.16774070067169,
          5.04632856830695,
          4.94726383191083,
          4.58770826515446,
          4.13835014251084,
          4.06045359796882,
          4.63530665635318,
          4.62485692334614,
          3.93341992672261,
          4.88762227733767,
          4.96825190546972,
          3.80832003906478,
          4.96057466313368,
          5.0939015418123,
          4.23984554358195,
          4.2685504137767,
          3.28500783834613,
          3.60201340344135,
          3.4002608223214,
          4.23362208987531,
          4.59490650795107,
          4.99824283163569,
          5.03285949305828
        ),
        Group = c(
          "A",
          "A",
          "A",
          "B",
          "B",
          "B",
          "C",
          "C",
          "C",
          "A",
          "A",
          "A",
          "B",
          "B",
          "B",
          "C",
          "C",
          "C",
          "A",
          "A",
          "A",
          "B",
          "B",
          "B",
          "C",
          "C",
          "C",
          "A",
          "A",
          "A",
          "B",
          "B",
          "B",
          "C",
          "C",
          "C",
          "A",
          "A",
          "A",
          "B",
          "B",
          "B",
          "C",
          "C",
          "C",
          "A",
          "A",
          "A",
          "B",
          "B",
          "B",
          "C",
          "C",
          "C"
        ),
        Replicate = c(
          1,
          2,
          3,
          1,
          2,
          3,
          1,
          2,
          3,
          1,
          2,
          3,
          1,
          2,
          3,
          1,
          2,
          3,
          1,
          2,
          3,
          1,
          2,
          3,
          1,
          2,
          3,
          1,
          2,
          3,
          1,
          2,
          3,
          1,
          2,
          3,
          1,
          2,
          3,
          1,
          2,
          3,
          1,
          2,
          3,
          1,
          2,
          3,
          1,
          2,
          3,
          1,
          2,
          3
        ),
        Batch = c(
          1,
          2,
          2,
          1,
          1,
          2,
          1,
          2,
          2,
          1,
          2,
          2,
          1,
          1,
          2,
          1,
          2,
          2,
          1,
          2,
          2,
          1,
          1,
          2,
          1,
          2,
          2,
          1,
          2,
          2,
          1,
          1,
          2,
          1,
          2,
          2,
          1,
          2,
          2,
          1,
          1,
          2,
          1,
          2,
          2,
          1,
          2,
          2,
          1,
          1,
          2,
          1,
          2,
          2
        ),
        Label = c(
          "A1",
          "A2",
          "A3",
          "B1",
          "B2",
          "B3",
          "C1",
          "C2",
          "C3",
          "A1",
          "A2",
          "A3",
          "B1",
          "B2",
          "B3",
          "C1",
          "C2",
          "C3",
          "A1",
          "A2",
          "A3",
          "B1",
          "B2",
          "B3",
          "C1",
          "C2",
          "C3",
          "A1",
          "A2",
          "A3",
          "B1",
          "B2",
          "B3",
          "C1",
          "C2",
          "C3",
          "A1",
          "A2",
          "A3",
          "B1",
          "B2",
          "B3",
          "C1",
          "C2",
          "C3",
          "A1",
          "A2",
          "A3",
          "B1",
          "B2",
          "B3",
          "C1",
          "C2",
          "C3"
        )
      ),
      row.names = c(NA, -54L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
})

test_that("plot_histogram resolves group colors by first observed group order", {
  color_values <- c("#5954d6", "#e1562c", "#b80058")
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
  plot <- plot_histogram(
    counts_dat,
    sample_metadata = nidap_sample_metadata,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    group_colname = "Group",
    color_values = color_values,
    color_by_group = TRUE
  )
  scales <- ggplot2::ggplot_build(plot)$plot$scales$scales
  colour_scale <- scales[[which(vapply(
    scales,
    function(scale) "colour" %in% scale$aesthetics,
    logical(1)
  ))[[1]]]]

  expect_equal(
    colour_scale$palette.cache,
    c(B = "#5954d6", A = "#e1562c", C = "#b80058")
  )
})

test_that("plot_histogram resolves group colors by factor level order", {
  color_values <- c("#5954d6", "#e1562c", "#b80058")
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

  plot <- plot_histogram(
    counts_dat,
    sample_metadata = sample_metadata,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    group_colname = "Group",
    color_values = color_values,
    color_by_group = TRUE
  )
  scales <- ggplot2::ggplot_build(plot)$plot$scales$scales
  colour_scale <- scales[[which(vapply(
    scales,
    function(scale) "colour" %in% scale$aesthetics,
    logical(1)
  ))[[1]]]]

  expect_equal(
    colour_scale$palette.cache,
    c(C = "#5954d6", A = "#e1562c", B = "#b80058")
  )
})

test_that("plot_histogram resolves sample colors by first observed sample order", {
  color_values <- c("#5954d6", "#e1562c", "#b80058", "#00c6f8")
  counts_dat <- nidap_filtered_counts[, c("Gene", "B1", "A1", "C1", "A2")]
  plot <- plot_histogram(
    counts_dat,
    sample_metadata = nidap_sample_metadata,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    color_values = color_values,
    color_by_group = FALSE
  )
  scales <- ggplot2::ggplot_build(plot)$plot$scales$scales
  colour_scale <- scales[[which(vapply(
    scales,
    function(scale) "colour" %in% scale$aesthetics,
    logical(1)
  ))[[1]]]]

  expect_equal(
    colour_scale$palette.cache,
    c(B1 = "#5954d6", A1 = "#e1562c", C1 = "#b80058", A2 = "#00c6f8")
  )
})

test_that("plot_histogram automatically sets log2 axis breaks when requested", {
  counts_dat <- data.frame(
    Gene = "gene_a",
    S1 = 0,
    S2 = 4,
    S3 = 16,
    check.names = FALSE
  )
  sample_metadata <- data.frame(
    Sample = c("S1", "S2", "S3"),
    Group = c("A", "A", "A"),
    check.names = FALSE
  )

  plot <- plot_histogram(
    counts_dat,
    sample_metadata = sample_metadata,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    color_by_group = FALSE,
    use_log2_x_axis = TRUE
  )
  x_scale <- plot$scales$get_scales("x")
  axis_breaks <- x_scale$breaks(c(0.5, 64.5))

  expect_equal(plot$data$count, c(0.5, 4.5, 16.5))
  expect_true(all(c(0.5, 4.5, 16.5, 64.5) %in% axis_breaks))
  expect_equal(x_scale$labels(c(0.5, 4.5, 16.5)), c("0", "4", "16"))
})

test_that("plot_histogram automatic log2 axis starts at zero", {
  counts_dat <- data.frame(
    Gene = "gene_a",
    S1 = 4,
    S2 = 16,
    check.names = FALSE
  )
  sample_metadata <- data.frame(
    Sample = c("S1", "S2"),
    Group = c("A", "A"),
    check.names = FALSE
  )

  plot <- plot_histogram(
    counts_dat,
    sample_metadata = sample_metadata,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    color_by_group = FALSE,
    use_log2_x_axis = TRUE
  )
  x_scale <- plot$scales$get_scales("x")

  expect_equal(2^x_scale$limits[1], 0.5)
  expect_equal(x_scale$labels(2^x_scale$limits[1]), "0")
})

test_that("plot_histogram works with tibbles", {
  p <- plot_histogram(
    nidap_filtered_counts,
    sample_metadata = nidap_sample_metadata,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
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
    color_by_group = FALSE,
    set_min_max_for_x_axis = FALSE,
    minimum_for_x_axis = -1,
    maximum_for_x_axis = 1,
    legend_position = "top",
    legend_font_size = 10,
    number_of_legend_columns = 6
  )
  expect_s3_class(p$layers[[1]], "ggproto")
  expect_s3_class(p$layers[[1]]$geom, "GeomArea")
  expect_equal(
    head(p$data),
    structure(
      list(
        Gene = c(
          "0610007P14Rik",
          "0610007P14Rik",
          "0610007P14Rik",
          "0610007P14Rik",
          "0610007P14Rik",
          "0610007P14Rik"
        ),
        Sample = c("A1", "A2", "A3", "B1", "B2", "B3"),
        count = c(1049, 950, 934, 1068, 1140, 947),
        Group = c("A", "A", "A", "B", "B", "B"),
        Replicate = c(1, 2, 3, 1, 2, 3),
        Batch = c(1, 2, 2, 1, 1, 2),
        Label = c("A1", "A2", "A3", "B1", "B2", "B3")
      ),
      row.names = c(NA, -6L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
})

test_that("plot_histogram result is the same for MOO or dataframe", {
  moo <- multiOmicDataSet(
    sample_metadata = nidap_sample_metadata,
    anno_dat = data.frame(),
    counts_lst = list("raw" = nidap_raw_counts)
  )
  expect_equal(
    plot_histogram(moo, count_type = "raw"),
    plot_histogram(nidap_raw_counts, sample_meta = nidap_sample_metadata)
  )
})

test_that("plot_histogram accepts print_plots and save_plots via moo dispatch without error", {
  moo <- multiOmicDataSet(
    sample_metadata = nidap_sample_metadata,
    anno_dat = data.frame(),
    counts_lst = list("raw" = nidap_raw_counts)
  )
  expect_no_error(
    plot_histogram(
      moo,
      count_type = "raw",
      group_colname = "Group",
      color_by_group = TRUE,
      print_plots = FALSE,
      save_plots = FALSE
    )
  )
})
