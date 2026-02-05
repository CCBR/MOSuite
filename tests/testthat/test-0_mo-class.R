test_that("constructing MOO works for RENEE data", {
  moo <- create_multiOmicDataSet_from_files(
    system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite"),
    system.file(
      "extdata",
      "RSEM.genes.expected_count.all_samples.txt.gz",
      package = "MOSuite"
    ),
    sample_id_colname = "sample_id",
    feature_id_colname = "gene_id"
  )
  expect_equal(
    moo@sample_meta,
    structure(
      list(
        sample_id = c("KO_S3", "KO_S4", "WT_S1", "WT_S2"),
        condition = c("knockout", "knockout", "wildtype", "wildtype")
      ),
      row.names = c(NA, -4L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
  expect_equal(
    moo@annotation |> head(),
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
        GeneName = c("A1BG", "A1BG-AS1", "A1CF", "A2M", "A2M-AS1", "A2ML1")
      ),
      row.names = c(NA, -6L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
  expect_equal(
    moo@counts$raw |> head(),
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
      row.names = c(NA, -6L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
})

test_that("constructing MOO works from CSV files", {
  moo <- create_multiOmicDataSet_from_files(
    system.file(
      "extdata",
      "nidap",
      "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz",
      package = "MOSuite"
    ),
    system.file("extdata", "nidap", "Raw_Counts.csv.gz", package = "MOSuite"),
    delim = ","
  )
  expect_equal(
    moo@sample_meta,
    structure(
      list(
        Sample = c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"),
        Group = c("A", "A", "A", "B", "B", "B", "C", "C", "C"),
        Replicate = c(1, 2, 3, 1, 2, 3, 1, 2, 3),
        Batch = c(1, 2, 2, 1, 1, 2, 1, 2, 2),
        Label = c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")
      ),
      row.names = c(NA, -9L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
  expect_equal(
    moo@annotation |> head(),
    structure(
      list(
        GeneName = c(
          "RP23-271O17.1",
          "Gm26206",
          "Xkr4",
          "RP23-317L18.1",
          "RP23-317L18.4",
          "RP23-317L18.3"
        )
      ),
      row.names = c(NA, -6L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
  expect_equal(
    moo@counts$raw |> head(),
    structure(
      list(
        GeneName = c(
          "RP23-271O17.1",
          "Gm26206",
          "Xkr4",
          "RP23-317L18.1",
          "RP23-317L18.4",
          "RP23-317L18.3"
        ),
        A1 = c(0, 0, 0, 0, 0, 0),
        A2 = c(0, 0, 0, 0, 0, 0),
        A3 = c(0, 0, 0, 0, 0, 0),
        B1 = c(0, 0, 0, 0, 0, 0),
        B2 = c(0, 0, 0, 0, 0, 0),
        B3 = c(0, 0, 0, 0, 0, 0),
        C1 = c(0, 0, 0, 0, 0, 0),
        C2 = c(0, 0, 0, 0, 0, 0),
        C3 = c(0, 0, 0, 0, 0, 0)
      ),
      row.names = c(NA, -6L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
})

test_that("annotation minimally contains feature id column", {
  moo <- create_multiOmicDataSet_from_dataframes(
    readr::read_tsv(
      system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite")
    ),
    gene_counts |> glue_gene_symbols()
  )
  expect_equal(
    moo@annotation |> head(),
    structure(
      list(
        gene_id = structure(
          c(
            "ENSG00000121410.11|A1BG",
            "ENSG00000268895.5|A1BG-AS1",
            "ENSG00000148584.15|A1CF",
            "ENSG00000175899.14|A2M",
            "ENSG00000245105.3|A2M-AS1",
            "ENSG00000166535.20|A2ML1"
          ),
          class = c("glue", "character")
        )
      ),
      row.names = c(NA, -6L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
})

test_that("multiOmicDataSet from data frames detect problems", {
  sample_meta <- data.frame(
    sample_id = c("KO_S3", "KO_S4", "WT_S1", "WT_S2"),
    condition = factor(
      c("knockout", "knockout", "wildtype", "wildtype"),
      levels = c("wildtype", "knockout")
    )
  )
  expect_error(
    create_multiOmicDataSet_from_dataframes(sample_meta, gene_counts[, 1:4]),
    "Not all sample IDs in the sample metadata are in the count data"
  )
})

test_that("extract_counts works", {
  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts),
      "norm" = list("voom" = as.data.frame(nidap_norm_counts))
    )
  )
  expect_equal(extract_counts(moo, "clean"), moo@counts$clean)
  expect_equal(extract_counts(moo, "norm", "voom"), moo@counts$norm$voom)
  expect_error(extract_counts(moo, "notacounttype"), "not in moo")
  expect_error(
    extract_counts(moo, "raw", "notasubtype"),
    "does not contain subtypes"
  )
  expect_error(extract_counts(moo, "norm"), "contains subtypes")
})


test_that("write_multiOmicDataSet_properties works", {
  moo_nidap <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts),
      "norm" = list("voom" = as.data.frame(nidap_norm_counts))
    )
  ) |>
    diff_counts(
      count_type = "filt",
      sub_count_type = NULL,
      sample_id_colname = "Sample",
      feature_id_colname = "Gene",
      covariates_colnames = c("Group", "Batch"),
      contrast_colname = c("Group"),
      contrasts = c("B-A", "C-A", "B-C"),
      voom_normalization_method = "quantile",
    ) |>
    filter_diff(
      significance_column = "adjpval",
      significance_cutoff = 0.05,
      change_column = "logFC",
      change_cutoff = 1,
      filtering_mode = "any",
      include_estimates = c("FC", "logFC", "tstat", "pval", "adjpval"),
      round_estimates = TRUE,
      rounding_decimal_for_percent_cells = 0,
      contrast_filter = "none",
      contrasts = c(),
      groups = c(),
      groups_filter = "none",
      label_font_size = 6,
      label_distance = 1,
      y_axis_expansion = 0.08,
      fill_colors = c("steelblue1", "whitesmoke"),
      pie_chart_in_3d = TRUE,
      bar_width = 0.4,
      draw_bar_border = TRUE,
      plot_type = "bar",
      plot_titles_fontsize = 12
    )
  moo_nidap@analyses$foo <- "bar"

  temp_dir <- tempfile(pattern = "moo-write-")
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  expect_equal(write_multiOmicDataSet_properties(moo_nidap, temp_dir), temp_dir)

  expect_true(file.exists(file.path(temp_dir, "sample_metadata.csv")))
  expect_true(file.exists(file.path(temp_dir, "feature_annotation.csv")))

  expect_true(file.exists(file.path(temp_dir, "counts", "raw_counts.csv")))
  expect_true(file.exists(file.path(temp_dir, "counts", "clean_counts.csv")))
  expect_true(file.exists(file.path(temp_dir, "counts", "filt_counts.csv")))
  expect_true(
    file.exists(file.path(temp_dir, "counts", "norm", "voom_counts.csv"))
  )

  expect_true(file.exists(file.path(temp_dir, "analyses", "foo.rds")))
  expect_true(file.exists(file.path(
    temp_dir,
    "analyses",
    "diff",
    "diff_B-A.csv"
  )))
  expect_true(file.exists(file.path(
    temp_dir,
    "analyses",
    "diff",
    "diff_C-A.csv"
  )))
  expect_true(file.exists(file.path(
    temp_dir,
    "analyses",
    "diff",
    "diff_B-C.csv"
  )))
  expect_true(file.exists(file.path(temp_dir, "analyses", "diff_filt.csv")))
  expect_true(dir.exists(file.path(temp_dir, "analyses", "colors")))
  expect_true(file.exists(file.path(
    temp_dir,
    "analyses",
    "colors",
    "colors_Sample.rds"
  )))
  expect_true(file.exists(file.path(
    temp_dir,
    "analyses",
    "colors",
    "colors_Batch.rds"
  )))
  expect_true(file.exists(file.path(
    temp_dir,
    "analyses",
    "colors",
    "colors_Group.rds"
  )))
  expect_true(file.exists(file.path(
    temp_dir,
    "analyses",
    "colors",
    "colors_Label.rds"
  )))
  expect_true(file.exists(file.path(
    temp_dir,
    "analyses",
    "colors",
    "colors_Replicate.rds"
  )))
})

test_that("write_multiOmicDataSet and read_multiOmicDataSet work", {
  # Create a simple moo object
  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(GeneName = unique(nidap_raw_counts$GeneName)),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts)
    )
  )

  # Write to temp file
  temp_file <- tempfile(pattern = "moo-", fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  # Test write returns filepath invisibly
  expect_equal(write_multiOmicDataSet(moo, temp_file), temp_file)
  expect_true(file.exists(temp_file))

  # Test read
  moo_read <- read_multiOmicDataSet(temp_file)
  expect_true(S7::S7_inherits(moo_read, multiOmicDataSet))

  # Verify all properties match
  expect_equal(moo_read@sample_meta, moo@sample_meta)
  expect_equal(moo_read@annotation, moo@annotation)
  expect_equal(moo_read@counts, moo@counts)
  expect_equal(names(moo_read@analyses), names(moo@analyses))
})

test_that("write_multiOmicDataSet validates input", {
  expect_error(
    write_multiOmicDataSet("not a moo"),
    "moo must be a multiOmicDataSet"
  )
  expect_error(
    write_multiOmicDataSet(list(sample_meta = data.frame())),
    "moo must be a multiOmicDataSet"
  )
})

test_that("read_multiOmicDataSet validates input", {
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  # Write a non-moo object
  readr::write_rds(list(a = 1, b = 2), temp_file)

  expect_error(
    read_multiOmicDataSet(temp_file),
    "RDS does not contain a multiOmicDataSet"
  )
})

test_that("write and read preserves complex moo with analyses", {
  moo_complex <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts),
      "norm" = list("voom" = as.data.frame(nidap_norm_counts))
    )
  )

  temp_file <- tempfile(pattern = "moo-complex-", fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)

  write_multiOmicDataSet(moo_complex, temp_file)
  moo_restored <- read_multiOmicDataSet(temp_file)

  expect_equal(moo_restored@sample_meta, moo_complex@sample_meta)
  expect_equal(moo_restored@annotation, moo_complex@annotation)
  expect_equal(moo_restored@counts$raw, moo_complex@counts$raw)
  expect_equal(moo_restored@counts$clean, moo_complex@counts$clean)
  expect_equal(moo_restored@counts$filt, moo_complex@counts$filt)
  expect_equal(moo_restored@counts$norm$voom, moo_complex@counts$norm$voom)
  expect_equal(
    names(moo_restored@analyses$colors),
    names(moo_complex@analyses$colors)
  )
})

test_that("validator returns character vector for invalid objects", {
  # Test validator with invalid count type
  expect_error(
    multiOmicDataSet(
      sample_metadata = as.data.frame(nidap_sample_metadata),
      anno_dat = data.frame(),
      counts_lst = list(
        "raw" = as.data.frame(nidap_raw_counts),
        "invalid_type" = as.data.frame(nidap_clean_raw_counts)
      )
    ),
    "@counts can only contain these names"
  )

  # Test validator with missing raw counts
  expect_error(
    multiOmicDataSet(
      sample_metadata = as.data.frame(nidap_sample_metadata),
      anno_dat = data.frame(),
      counts_lst = list(
        "clean" = as.data.frame(nidap_clean_raw_counts)
      )
    ),
    "@counts must contain at least 'raw' counts"
  )

  # Test validator with mismatched sample IDs
  mismatched_counts <- as.data.frame(nidap_raw_counts)
  colnames(mismatched_counts)[2] <- "WRONG_SAMPLE_ID"
  expect_error(
    multiOmicDataSet(
      sample_metadata = as.data.frame(nidap_sample_metadata),
      anno_dat = data.frame(),
      counts_lst = list(
        "raw" = mismatched_counts
      )
    ),
    "@sample_meta"
  )
})

test_that("validator returns NULL for valid objects", {
  # Create a valid moo object
  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts)
    )
  )
  # If the validator returned errors, the object wouldn't have been created
  # So we just check that the object exists and is the correct type
  expect_true(S7::S7_inherits(moo, multiOmicDataSet))
})
