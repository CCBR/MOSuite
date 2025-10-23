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
  sample_meta <- system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite") %>%
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
  moo2 <- moo %>% set_color_pal(
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
