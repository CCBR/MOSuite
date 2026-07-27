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
        `1` = "#00c6f8",
        `2` = "#d163e6",
        `3` = "#00a76c"
      ),
      Batch = c(`1` = "#ff9287", `2` = "#008cf9"),
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
test_that("get_colors_lst handles alternative palette vectors", {
  sample_meta <- system.file(
    "extdata",
    "sample_metadata.tsv.gz",
    package = "MOSuite"
  ) |>
    readr::read_tsv()
  result <- get_colors_lst(
    sample_meta,
    palette = RColorBrewer::brewer.pal(12, "Set3")
  )
  expect_type(result, "list")
  expect_length(result, ncol(sample_meta))
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

test_that("get_colors_vctr retries from offset 0 when offset pushes past palette end", {
  # Palette has 12 colors; offset=10, n_obs=5: 10+5=15 > 12, so retry without offset.
  # The column should receive palette colors 1-5, not random colors.
  dat <- data.frame(group = paste0("x", seq_len(5)))
  result <- get_colors_vctr(dat, "group", color_offset = 10L)
  expect_length(result, 5)
  expect_named(result, paste0("x", seq_len(5)))
  # Should get the first 5 palette colors (same as offset=0), not random
  expected <- get_colors_vctr(dat, "group", color_offset = 0L)
  expect_equal(result, expected)
})

test_that("get_colors_lst columns exceeding palette size fall back to random colors", {
  # 13 unique values exceeds the 12-color mosuite_palette; should message and use random colors
  dat_big <- data.frame(group = paste0("cat", seq_len(13)))
  expect_no_warning(
    expect_message(
      result <- get_colors_lst(dat_big),
      "exceeds the palette maximum"
    )
  )
  expect_length(result$group, 13)
  expect_named(result$group, paste0("cat", seq_len(13)))
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

test_that("select_mosuite_colors returns n colors from mosuite_palette", {
  result <- select_mosuite_colors(3)
  expect_length(result, 3)
  expect_true(all(grepl("^#", result)))
})

test_that("select_mosuite_colors clamps to palette length when n exceeds it", {
  pal_len <- length(mosuite_palette)
  result <- select_mosuite_colors(pal_len + 10)
  expect_length(result, pal_len)
})

test_that("get_observed_values drops NAs and returns unique values in first-seen order", {
  dat <- data.frame(group = c("B", NA, "A", "B", NA))
  result <- MOSuite:::get_observed_values(dat, "group")
  expect_equal(result, c("B", "A"))
})

test_that("get_observed_values respects factor level order and excludes unseen levels", {
  dat <- data.frame(
    group = factor(c("B", "A", NA), levels = c("C", "A", "B", "D"))
  )
  result <- MOSuite:::get_observed_values(dat, "group")
  # "C" and "D" are levels but never observed; should be excluded
  expect_equal(result, c("A", "B"))
})

test_that("get_colors_vctr returns empty character vector when column has zero observed values", {
  expect_equal(
    get_colors_vctr(data.frame(group = character(0)), "group"),
    character(0)
  )
  expect_equal(
    get_colors_vctr(
      data.frame(group = c(NA_character_, NA_character_)),
      "group"
    ),
    character(0)
  )
})

test_that("get_colors_vctr handles a column with NA values", {
  dat <- data.frame(group = c("A", NA, "B", "A"))
  result <- get_colors_vctr(dat, "group")
  expect_length(result, 2)
  expect_named(result, c("A", "B"))
})

test_that("resolve_plot_colors returns color_values unchanged when column has no observations", {
  dat <- data.frame(group = character(0))
  colors <- c(A = "red")
  result <- MOSuite:::resolve_plot_colors(dat, "group", color_values = colors)
  expect_equal(result, colors)
})

test_that("resolve_plot_colors returns color_values unchanged when column is all NA", {
  dat <- data.frame(group = c(NA_character_, NA_character_))
  colors <- c(A = "red")
  result <- MOSuite:::resolve_plot_colors(dat, "group", color_values = colors)
  expect_equal(result, colors)
})

test_that("display_palette returns a ggplot invisibly", {
  result <- display_palette(c("#FF0000", "#00FF00", "#0000FF"))
  expect_s3_class(result, "gg")
})

test_that("display_colors returns a patchwork object", {
  moo <- create_multiOmicDataSet_from_dataframes(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    counts_dat = as.data.frame(nidap_raw_counts)
  )
  result <- display_colors(moo)
  expect_s3_class(result, "patchwork")
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
        `1` = "#00c6f8",
        `2` = "#d163e6",
        `3` = "#00a76c"
      ),
      Batch = c(`1` = "#ff9287", `2` = "#008cf9"),
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
      palette = RColorBrewer::brewer.pal(3, "Set2")
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
        `1` = "#00c6f8",
        `2` = "#d163e6",
        `3` = "#00a76c"
      ),
      Batch = c(`1` = "#ff9287", `2` = "#008cf9"),
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
