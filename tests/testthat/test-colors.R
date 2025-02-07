test_that("get_random_colors works", {
  set.seed(10)
  expect_equal(
    get_random_colors(5),
    c("#B85CD0", "#B4E16D", "#DC967D", "#A6DCC5", "#B5AAD3")
  )
  expect_equal(
    get_random_colors(3),
    c("#B3C4C7", "#B7D579", "#C56BC8")
  )
  expect_error(get_random_colors(0), "num_colors must be at least 1")
})

test_that("set_colors works on nidap_sample_metadata", {
  expect_equal(
    set_colors(nidap_sample_metadata),
    list(
      Sample = c(
        A1 = "#000000", A2 = "#E69F00", A3 = "#56B4E9",
        B1 = "#009E73", B2 = "#F0E442", B3 = "#0072B2", C1 = "#D55E00",
        C2 = "#CC79A7", C3 = "#999999"
      ), Group = c(
        A = "#000000", B = "#E69F00",
        C = "#56B4E9"
      ), Replicate = c(
        `1` = "#000000", `2` = "#E69F00",
        `3` = "#56B4E9"
      ), Batch = c(`1` = "#000000", `2` = "#E69F00"),
      Label = c(
        A1 = "#000000", A2 = "#E69F00", A3 = "#56B4E9",
        B1 = "#009E73", B2 = "#F0E442", B3 = "#0072B2", C1 = "#D55E00",
        C2 = "#CC79A7", C3 = "#999999"
      )
    )
  )
})
test_that("set_colors handles alternative palette functions", {
  sample_meta <- system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite") %>%
    readr::read_tsv()
  expect_message(
    expect_warning(
      set_colors(sample_meta, palette_fun = RColorBrewer::brewer.pal, name = "Set3"),
      "minimal value for n is 3"
    ),
    "Warning raised in set_colors"
  )
})
