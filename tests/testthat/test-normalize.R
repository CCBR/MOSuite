test_that("normalize works for NIDAP", {
  moo <- multiOmicDataSet(
    sample_meta_dat = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts)
    )
  ) %>%
    normalize_counts(group_column = "Group", label_column = "Label")
  expect_true(equal_dfs(
    moo@counts[["norm"]][["voom"]] %>%
      dplyr::arrange(desc(Gene)),
    as.data.frame(nidap_norm_counts) %>%
      dplyr::arrange(desc(Gene))
  ))
})

test_that("noramlize works for RENEE", {
  moo <- create_multiOmicDataSet_from_dataframes(readr::read_tsv(
    system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite")
  ), gene_counts) %>%
    clean_raw_counts() %>%
    filter_counts(
      group_column = "condition",
      label_column = "sample_id",
      minimum_count_value_to_be_considered_nonzero = 1,
      minimum_number_of_samples_with_nonzero_counts_in_total = 1,
      minimum_number_of_samples_with_nonzero_counts_in_a_group = 1,
      make_plots = FALSE
    ) %>%
    normalize_counts(group_column = "condition", label_column = "sample_id")
  expect_equal(
    head(moo@counts$norm$voom),
    structure(
      list(
        gene_id = c(
          "ENSG00000215458.8",
          "ENSG00000160179.18",
          "ENSG00000258017.1",
          "ENSG00000282393.1",
          "ENSG00000286104.1",
          "ENSG00000274422.1"
        ),
        KO_S3 = c(
          11.0751960068561,
          9.6086338540783,
          9.6086338540783,
          8.81615260371772,
          9.6086338540783,
          8.81615260371772
        ),
        KO_S4 = c(
          12.3480907442867,
          12.7703165561761,
          8.81615260371772,
          9.6086338540783,
          8.81615260371772,
          9.6086338540783
        ),
        WT_S1 = c(
          8.81615260371772,
          12.3480907442867,
          8.81615260371772,
          8.81615260371772,
          8.81615260371772,
          8.81615260371772
        ),
        WT_S2 = c(
          10.0048744792586,
          12.2369960496953,
          8.81615260371772,
          8.81615260371772,
          8.81615260371772,
          8.81615260371772
        )
      ),
      row.names = c(NA, 6L),
      class = "data.frame"
    )
  )
  expect_equal(
    tail(moo@counts$norm$voom),
    structure(
      list(
        gene_id = c(
          "ENSG00000157538.14",
          "ENSG00000160193.11",
          "ENSG00000182093.15",
          "ENSG00000182362.14",
          "ENSG00000173276.14",
          "ENSG00000237232.7"
        ),
        KO_S3 = c(
          12.3480907442867,
          9.6086338540783,
          11.8597009422769,
          11.0751960068561,
          11.8597009422769,
          8.81615260371772
        ),
        KO_S4 = c(
          12.7703165561761,
          9.6086338540783,
          9.6086338540783,
          8.81615260371772,
          12.7703165561761,
          9.6086338540783
        ),
        WT_S1 = c(
          12.2426956580003,
          10.5853565029804,
          11.7865266999202,
          8.81615260371772,
          11.7865266999202,
          8.81615260371772
        ),
        WT_S2 = c(
          12.4720029602325,
          10.9249210977479,
          11.4357186116131,
          8.81615260371772,
          12.2369960496953,
          8.81615260371772
        )
      ),
      row.names = 286:291,
      class = "data.frame"
    )
  )
})
