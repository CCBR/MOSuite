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
  expect_equal(moo@sample_meta, structure(
    list(
      sample_id = c("KO_S3", "KO_S4", "WT_S1", "WT_S2"),
      condition = c("knockout", "knockout", "wildtype", "wildtype")
    ),
    row.names = c(NA, -4L),
    class = c("tbl_df", "tbl", "data.frame")
  ))
  expect_equal(
    moo@annotation %>% head(),
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
    moo@counts$raw %>% head(),
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
  moo <- create_multiOmicDataSet_from_files(system.file("extdata", "nidap", "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz", package = "MOSuite"),
    system.file("extdata", "nidap", "Raw_Counts.csv.gz", package = "MOSuite"),
    delim = ","
  )
  expect_equal(moo@sample_meta, structure(
    list(
      Sample = c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"),
      Group = c("A", "A", "A", "B", "B", "B", "C", "C", "C"),
      Replicate = c(1, 2, 3, 1, 2, 3, 1, 2, 3),
      Batch = c(1, 2, 2, 1, 1, 2, 1, 2, 2),
      Label = c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")
    ),
    row.names = c(NA, -9L),
    class = c("tbl_df", "tbl", "data.frame")
  ))
  expect_equal(
    moo@annotation %>% head(),
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
    moo@counts$raw %>% head(),
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
    gene_counts %>% glue_gene_symbols()
  )
  expect_equal(
    moo@annotation %>% head(),
    structure(
      list(gene_id = structure(
        c(
          "ENSG00000121410.11|A1BG",
          "ENSG00000268895.5|A1BG-AS1",
          "ENSG00000148584.15|A1CF",
          "ENSG00000175899.14|A2M",
          "ENSG00000245105.3|A2M-AS1",
          "ENSG00000166535.20|A2ML1"
        ),
        class = c("glue", "character")
      )),
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
      "norm" = list(
        "voom" = as.data.frame(nidap_norm_counts)
      )
    )
  )
  expect_equal(extract_counts(moo, "clean"), moo@counts$clean)
  expect_equal(extract_counts(moo, "norm", "voom"), moo@counts$norm$voom)
  expect_error(extract_counts(moo, "notacounttype"), "not in moo")
  expect_error(extract_counts(moo, "raw", "notasubtype"), "does not contain subtypes")
  expect_error(extract_counts(moo, "norm"), "contains subtypes")
})
