test_that("get_random_colors works", {
  set.seed(10)
  expect_equal(
    get_random_colors(5),
    c("#B85CD0", "#B4E16D", "#DC967D", "#A6DCC5", "#B5AAD3")
  )
  expect_equal(get_random_colors(3), c("#B3C4C7", "#B7D579", "#C56BC8"))
  expect_error(get_random_colors(0), "num_colors must be at least 1")
})

test_that("get_colors_lst works on nidap_sample_metadata", {
  expect_equal(
    get_colors_lst(nidap_sample_metadata),
    list(
      Sample = c(
        A1 = "#000000",
        A2 = "#E69F00",
        A3 = "#56B4E9",
        B1 = "#009E73",
        B2 = "#F0E442",
        B3 = "#0072B2",
        C1 = "#D55E00",
        C2 = "#CC79A7",
        C3 = "#999999"
      ),
      Group = c(
        A = "#000000",
        B = "#E69F00",
        C = "#56B4E9"
      ),
      Replicate = c(
        `1` = "#000000",
        `2` = "#E69F00",
        `3` = "#56B4E9"
      ),
      Batch = c(`1` = "#000000", `2` = "#E69F00"),
      Label = c(
        A1 = "#000000",
        A2 = "#E69F00",
        A3 = "#56B4E9",
        B1 = "#009E73",
        B2 = "#F0E442",
        B3 = "#0072B2",
        C1 = "#D55E00",
        C2 = "#CC79A7",
        C3 = "#999999"
      )
    )
  )
})
test_that("get_colors_lst handles alternative palette functions", {
  sample_meta <- system.file(
    "extdata",
    "sample_metadata.tsv.gz",
    package = "MOSuite"
  ) |>
    readr::read_tsv()
  expect_message(
    expect_warning(
      get_colors_lst(
        sample_meta,
        palette_fun = RColorBrewer::brewer.pal,
        name = "Set3"
      ),
      "minimal value for n is 3"
    ),
    "Warning raised in "
  )
})
test_that("get_colors_vctr falls back to random colors when n exceeds palette max", {
  # Okabe-Ito palette has a maximum of 9 colors. When n > 9, the function
  # should fall back to get_random_colors() and emit a message.
  dat_many_cats <- data.frame(
    group = paste0("cat", seq_len(12))
  )
  expect_no_warning(
    expect_message(
      result <- get_colors_vctr(dat_many_cats, "group"),
      "exceeds the palette maximum"
    )
  )
  expect_length(result, 12)
  expect_named(result, paste0("cat", seq_len(12)))
})

test_that("resolve_plot_colors preserves named color mappings", {
  dat <- data.frame(group = c("B", "A", "C", "A"))
  colors <- c(A = "red", B = "blue", C = "green")

  expect_equal(resolve_plot_colors(dat, "group", colors), colors)
})

test_that("resolve_plot_colors names palettes by first observed category order", {
  dat <- data.frame(group = c("B", "A", "C", "A"))
  colors <- c("red", "blue", "green")

  expect_equal(
    resolve_plot_colors(dat, "group", colors),
    c(B = "red", A = "blue", C = "green")
  )
})

test_that("resolve_plot_colors generates colors when none are supplied", {
  dat <- data.frame(group = c("B", "A", "C", "A"))

  expect_equal(
    resolve_plot_colors(dat, "group"),
    c(B = "#000000", A = "#E69F00", C = "#56B4E9")
  )
})

test_that("resolve_plot_colors rejects too few explicit colors", {
  dat <- data.frame(group = c("B", "A", "C", "A"))

  expect_error(
    resolve_plot_colors(dat, "group", c("red", "blue")),
    "must contain at least 3 colors"
  )
})

test_that("resolve_plot_colors treats non-matching names as palette labels", {
  dat <- data.frame(group = c("B", "A", "C", "A"))

  expect_equal(
    resolve_plot_colors(dat, "group", c(indigo = "red", carrot = "blue", jade = "green")),
    c(B = "red", A = "blue", C = "green")
  )
})

test_that("set_color_pal overrides the color palette", {
  moo <- create_multiOmicDataSet_from_dataframes(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    counts_dat = as.data.frame(nidap_raw_counts)
  )
  expect_equal(
    moo@analyses$colors,
    list(
      Sample = c(
        A1 = "#000000",
        A2 = "#E69F00",
        A3 = "#56B4E9",
        B1 = "#009E73",
        B2 = "#F0E442",
        B3 = "#0072B2",
        C1 = "#D55E00",
        C2 = "#CC79A7",
        C3 = "#999999"
      ),
      Group = c(
        A = "#000000",
        B = "#E69F00",
        C = "#56B4E9"
      ),
      Replicate = c(
        `1` = "#000000",
        `2` = "#E69F00",
        `3` = "#56B4E9"
      ),
      Batch = c(`1` = "#000000", `2` = "#E69F00"),
      Label = c(
        A1 = "#000000",
        A2 = "#E69F00",
        A3 = "#56B4E9",
        B1 = "#009E73",
        B2 = "#F0E442",
        B3 = "#0072B2",
        C1 = "#D55E00",
        C2 = "#CC79A7",
        C3 = "#999999"
      )
    )
  )
  moo2 <- moo |>
    set_color_pal(
      colname = "Group",
      palette_fun = RColorBrewer::brewer.pal,
      name = "Set2"
    )
  expect_equal(
    moo2@analyses$colors,
    list(
      Sample = c(
        A1 = "#000000",
        A2 = "#E69F00",
        A3 = "#56B4E9",
        B1 = "#009E73",
        B2 = "#F0E442",
        B3 = "#0072B2",
        C1 = "#D55E00",
        C2 = "#CC79A7",
        C3 = "#999999"
      ),
      Group = c(
        A = "#66C2A5",
        B = "#FC8D62",
        C = "#8DA0CB"
      ),
      Replicate = c(
        `1` = "#000000",
        `2` = "#E69F00",
        `3` = "#56B4E9"
      ),
      Batch = c(`1` = "#000000", `2` = "#E69F00"),
      Label = c(
        A1 = "#000000",
        A2 = "#E69F00",
        A3 = "#56B4E9",
        B1 = "#009E73",
        B2 = "#F0E442",
        B3 = "#0072B2",
        C1 = "#D55E00",
        C2 = "#CC79A7",
        C3 = "#999999"
      )
    )
  )
})
