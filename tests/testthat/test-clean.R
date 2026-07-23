test_that("clean_raw_counts works for NIDAP data", {
  moo_nidap <- create_multiOmicDataSet_from_dataframes(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    counts_dat = as.data.frame(nidap_raw_counts)
  ) |>
    clean_raw_counts(print_plots = TRUE)

  actual <- moo_nidap@counts[["clean"]] |>
    dplyr::rename(Gene = GeneName) |>
    as.data.frame()

  expected <- as.data.frame(nidap_clean_raw_counts)

  cmp <- all.equal(actual, expected, check.attributes = FALSE)
  expect_true(isTRUE(cmp), info = paste(cmp, collapse = "\n"))
})

test_that("clean_raw_counts works for RENEE data", {
  moo <- create_multiOmicDataSet_from_dataframes(
    readr::read_tsv(
      system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite")
    ),
    gene_counts
  ) |>
    clean_raw_counts()
  expect_equal(
    head(moo@counts$clean),
    structure(
      list(
        gene_id = c(
          "ENSG00000121410.11",
          "ENSG00000268895.5",
          "ENSG00000148584.15",
          "ENSG00000175899.14",
          "ENSG00000245105.3",
          "ENSG00000166535.20"
        ),
        KO_S3 = c(0, 0, 0, 0, 0, 0),
        KO_S4 = c(0, 0, 0, 0, 0, 0),
        WT_S1 = c(0, 0, 0, 0, 0, 0),
        WT_S2 = c(0, 0, 0, 0, 0, 0)
      ),
      row.names = c(NA, 6L),
      class = "data.frame"
    )
  )
  expect_equal(
    tail(moo@counts$clean),
    structure(
      list(
        gene_id = c(
          "ENSG00000232242.2",
          "ENSG00000162378.13",
          "ENSG00000159840.16",
          "ENSG00000274572.1",
          "ENSG00000074755.15",
          "ENSG00000272920.1"
        ),
        KO_S3 = c(0, 0, 0, 0, 0, 0),
        KO_S4 = c(0, 0, 0, 0, 0, 0),
        WT_S1 = c(0, 0, 0, 0, 0, 0),
        WT_S2 = c(0, 0, 0, 0, 0, 0)
      ),
      row.names = 58924:58929,
      class = "data.frame"
    )
  )
})

test_that("clean_raw_counts plots CPM histogram after cleaning without all-zero genes", {
  captured_histogram_counts <- NULL
  histogram_args <- NULL
  saved_filenames <- character()
  group_colors <- c(A = "red", B = "blue")

  local_mocked_bindings(
    plot_histogram = function(moo_counts, ...) {
      captured_histogram_counts <<- moo_counts
      histogram_args <<- list(...)
      ggplot2::ggplot()
    },
    print_or_save_plot = function(plot, filename, ...) {
      saved_filenames <<- c(saved_filenames, basename(filename))
      invisible(NULL)
    },
    .package = "MOSuite"
  )

  counts_dat <- data.frame(
    GeneName = c("all_zero", "keep_one", "keep_one"),
    S1 = c(0, 10, 5),
    S2 = c(0, 0, 5),
    check.names = FALSE
  )
  sample_metadata <- data.frame(
    Sample = c("S1", "S2"),
    Group = c("A", "B"),
    check.names = FALSE
  )

  create_multiOmicDataSet_from_dataframes(
    sample_metadata = sample_metadata,
    counts_dat = counts_dat,
    sample_id_colname = "Sample",
    feature_id_colname = "GeneName"
  ) |>
    clean_raw_counts(
      sample_id_colname = "Sample",
      feature_id_colname = "GeneName",
      group_colname = "Group",
      colors_for_plots = group_colors,
      split_gene_name = FALSE,
      print_plots = TRUE,
      save_plots = FALSE
    )

  expect_true("read_depth.png" %in% saved_filenames)
  expect_true("cpm_histogram.png" %in% saved_filenames)
  expect_equal(captured_histogram_counts$GeneName, "keep_one")
  expect_equal(captured_histogram_counts$S1, 1e6)
  expect_equal(captured_histogram_counts$S2, 1e6)
  expect_equal(histogram_args$group_colname, "Group")
  expect_equal(histogram_args$color_values, group_colors)
  expect_true(histogram_args$color_by_group)
})

test_that("aggregate_duplicate_gene_names returns collapsed dfout", {
  counts_dat <- data.frame(
    gene_id = c("A", "A", "B"),
    sample1 = c(1, 2, 3),
    sample2 = c(4, 5, 6),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Case 1: aggregation enabled
  out <- MOSuite:::aggregate_duplicate_gene_names(
    counts_dat = counts_dat,
    gene_name_column_to_use_for_collapsing_duplicates = "gene_id",
    aggregate_rows_with_duplicate_gene_names = TRUE,
    split_gene_name = FALSE
  )

  expect_equal(nrow(out), 2)
  expect_equal(sum(duplicated(out$gene_id)), 0)

  a_row <- out[out$gene_id == "A", , drop = FALSE]
  expect_equal(a_row$sample1, 3)
  expect_equal(a_row$sample2, 9)

  # Case 2: aggregation disabled
  out_noagg <- MOSuite:::aggregate_duplicate_gene_names(
    counts_dat = counts_dat,
    gene_name_column_to_use_for_collapsing_duplicates = "gene_id",
    aggregate_rows_with_duplicate_gene_names = FALSE,
    split_gene_name = FALSE
  )

  expect_equal(nrow(out_noagg), 3)
  expect_equal(sum(duplicated(out_noagg$gene_id)), 1)
  expect_equal(out_noagg$sample1, counts_dat$sample1)
  expect_equal(out_noagg$sample2, counts_dat$sample2)
})
