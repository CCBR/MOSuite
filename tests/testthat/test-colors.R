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
        A1 = "#5954d6",
        A2 = "#e1562c",
        A3 = "#b80058",
        B1 = "#00c6f8",
        B2 = "#d163e6",
        B3 = "#00a76c",
        C1 = "#ff9287",
        C2 = "#008cf9",
        C3 = "#006e00"
      ),
      Group = c(
        A = "#5954d6",
        B = "#e1562c",
        C = "#b80058"
      ),
      Replicate = c(
        `1` = "#5954d6",
        `2` = "#e1562c",
        `3` = "#b80058"
      ),
      Batch = c(`1` = "#5954d6", `2` = "#e1562c"),
      Label = c(
        A1 = "#5954d6",
        A2 = "#e1562c",
        A3 = "#b80058",
        B1 = "#00c6f8",
        B2 = "#d163e6",
        B3 = "#00a76c",
        C1 = "#ff9287",
        C2 = "#008cf9",
        C3 = "#006e00"
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
  # MOSuite's default palette has 12 colors. When n > 12, the function
  # should fall back to get_random_colors() and emit a message.
  dat_many_cats <- data.frame(
    group = paste0("cat", seq_len(13))
  )
  expect_no_warning(
    expect_message(
      result <- get_colors_vctr(dat_many_cats, "group"),
      "exceeds the palette maximum"
    )
  )
  expect_length(result, 13)
  expect_named(result, paste0("cat", seq_len(13)))
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

test_that("color vectors use factor level order when grouping column is a factor", {
  dat <- data.frame(
    group = factor(c("B", "A", "C", "A"), levels = c("C", "A", "B", "D"))
  )

  expect_equal(
    get_colors_vctr(dat, "group"),
    c(C = "#5954d6", A = "#e1562c", B = "#b80058")
  )
  expect_equal(
    resolve_plot_colors(dat, "group", c("red", "blue", "green")),
    c(C = "red", A = "blue", B = "green")
  )
})

test_that("resolve_plot_colors generates colors when none are supplied", {
  dat <- data.frame(group = c("B", "A", "C", "A"))

  expect_equal(
    resolve_plot_colors(dat, "group"),
    c(B = "#5954d6", A = "#e1562c", C = "#b80058")
  )
})

test_that("resolve_plot_colors generates additional colors for too few explicit colors", {
  dat <- data.frame(group = c("B", "A", "C", "A"))

  expect_message(
    result <- resolve_plot_colors(dat, "group", c("red", "blue")),
    "Generating 1 additional colors"
  )
  expect_named(result, c("B", "A", "C"))
  expect_equal(unname(result[1:2]), c("red", "blue"))
  expect_equal(unname(result[3]), "#b80058")
})

test_that("resolve_plot_colors uses random fallback only through get_colors_vctr", {
  dat <- data.frame(group = paste0("cat", seq_len(13)))
  colors <- c(
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

  expect_message(
    expect_message(
      result <- resolve_plot_colors(dat, "group", colors),
      "Generating 1 additional colors"
    ),
    "exceeds the palette maximum"
  )
  expect_named(result, paste0("cat", seq_len(13)))
  expect_equal(unname(result[seq_along(colors)]), colors)
  expect_match(unname(result[13]), "^#[0-9A-F]{6}$")
})

test_that("resolve_plot_colors treats non-matching names as palette labels", {
  dat <- data.frame(group = c("B", "A", "C", "A"))

  expect_equal(
    resolve_plot_colors(
      dat,
      "group",
      c(indigo = "red", carrot = "blue", jade = "green")
    ),
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
        A1 = "#5954d6",
        A2 = "#e1562c",
        A3 = "#b80058",
        B1 = "#00c6f8",
        B2 = "#d163e6",
        B3 = "#00a76c",
        C1 = "#ff9287",
        C2 = "#008cf9",
        C3 = "#006e00"
      ),
      Group = c(
        A = "#5954d6",
        B = "#e1562c",
        C = "#b80058"
      ),
      Replicate = c(
        `1` = "#5954d6",
        `2` = "#e1562c",
        `3` = "#b80058"
      ),
      Batch = c(`1` = "#5954d6", `2` = "#e1562c"),
      Label = c(
        A1 = "#5954d6",
        A2 = "#e1562c",
        A3 = "#b80058",
        B1 = "#00c6f8",
        B2 = "#d163e6",
        B3 = "#00a76c",
        C1 = "#ff9287",
        C2 = "#008cf9",
        C3 = "#006e00"
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
        A1 = "#5954d6",
        A2 = "#e1562c",
        A3 = "#b80058",
        B1 = "#00c6f8",
        B2 = "#d163e6",
        B3 = "#00a76c",
        C1 = "#ff9287",
        C2 = "#008cf9",
        C3 = "#006e00"
      ),
      Group = c(
        A = "#66C2A5",
        B = "#FC8D62",
        C = "#8DA0CB"
      ),
      Replicate = c(
        `1` = "#5954d6",
        `2` = "#e1562c",
        `3` = "#b80058"
      ),
      Batch = c(`1` = "#5954d6", `2` = "#e1562c"),
      Label = c(
        A1 = "#5954d6",
        A2 = "#e1562c",
        A3 = "#b80058",
        B1 = "#00c6f8",
        B2 = "#d163e6",
        B3 = "#00a76c",
        C1 = "#ff9287",
        C2 = "#008cf9",
        C3 = "#006e00"
      )
    )
  )
})
