test_that("counts_dat_to_matrix works", {
  expect_equal(
    gene_counts |>
      dplyr::select(-GeneName) |>
      head() |>
      counts_dat_to_matrix(),
    structure(
      c(
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L,
        0L
      ),
      dim = c(6L, 4L),
      dimnames = list(
        c(
          "ENSG00000121410.11",
          "ENSG00000268895.5",
          "ENSG00000148584.15",
          "ENSG00000175899.14",
          "ENSG00000245105.3",
          "ENSG00000166535.20"
        ),
        c("KO_S3", "KO_S4", "WT_S1", "WT_S2")
      )
    )
  )
})

test_that("calc_cpm works on RENEE data", {
  sample_meta <- data.frame(
    sample_id = c("KO_S3", "KO_S4", "WT_S1", "WT_S2"),
    condition = factor(
      c("knockout", "knockout", "wildtype", "wildtype"),
      levels = c("wildtype", "knockout")
    )
  )
  moo <- create_multiOmicDataSet_from_dataframes(
    sample_meta,
    gene_counts |> dplyr::select(-GeneName)
  )
  moo <- moo |> calc_cpm()
  cpm_edger <- gene_counts |>
    dplyr::select(-GeneName) |>
    counts_dat_to_matrix() |>
    edgeR::cpm() |>
    as.data.frame() |>
    tibble::rownames_to_column("gene_id")
  expect_equal(moo@counts$cpm, cpm_edger)
})

test_that("calc_cpm_df works on NIDAP data", {
  df <- nidap_clean_raw_counts |> as.data.frame()
  trans.df <- df
  trans.df[, -1] <- edgeR::cpm(as.matrix(df[, -1]))

  expect_equal(
    calc_cpm_df(df, feature_id_colname = "Gene"),
    trans.df,
    ignore_attr = TRUE
  )
})
test_that("calc_cpm_df preserves rownames", {
  df <- nidap_clean_raw_counts |>
    as.data.frame() |>
    tail()
  trans.df <- df
  trans.df[, -1] <- edgeR::cpm(as.matrix(df[, -1]))

  expect_equal(
    calc_cpm_df(df, feature_id_colname = "Gene"),
    trans.df,
    ignore_attr = TRUE
  )
})

test_that("log_transform_counts transforms counts correctly", {
  dat <- data.frame(
    gene_id = c("g1", "g2"),
    s1 = c(1, 4),
    s2 = c(9, 0)
  )
  result <- log_transform_counts(dat, pseudocount = 1, base = 2)
  expect_equal(result$s1, log(c(1, 4) + 1, base = 2))
  expect_equal(result$s2, log(c(9, 0) + 1, base = 2))
})

test_that("log_transform_counts uses natural log by default", {
  dat <- data.frame(gene_id = "g1", s1 = 10)
  result <- log_transform_counts(dat, pseudocount = 0.5, base = "ln")
  expect_equal(result$s1, log(10.5))
})

test_that("log_transform_counts errors for invalid pseudocount", {
  dat <- data.frame(gene_id = "g1", s1 = 1)
  expect_error(
    log_transform_counts(dat, pseudocount = "a"),
    "pseudocount must be a single numeric value"
  )
  expect_error(
    log_transform_counts(dat, pseudocount = -1),
    "pseudocount cannot be negative"
  )
})

test_that("log_transform_counts errors when count + pseudocount <= 0", {
  dat <- data.frame(gene_id = "g1", s1 = -2)
  expect_error(log_transform_counts(dat, pseudocount = 1), "greater than 0")
})

test_that("resolve_log_transform_base returns exp(1) for string aliases", {
  expect_equal(resolve_log_transform_base("e"), exp(1))
  expect_equal(resolve_log_transform_base("ln"), exp(1))
  expect_equal(resolve_log_transform_base("natural"), exp(1))
  expect_equal(resolve_log_transform_base("E"), exp(1))
  expect_equal(resolve_log_transform_base("LN"), exp(1))
})

test_that("resolve_log_transform_base returns numeric base as-is", {
  expect_equal(resolve_log_transform_base(2), 2)
  expect_equal(resolve_log_transform_base(10), 10)
})

test_that("resolve_log_transform_base errors for invalid string", {
  expect_error(
    resolve_log_transform_base("log2"),
    "must be a single numeric value"
  )
})

test_that("resolve_log_transform_base errors for non-numeric non-string", {
  expect_error(
    resolve_log_transform_base(TRUE),
    "must be a single numeric value"
  )
  expect_error(
    resolve_log_transform_base(c(2, 10)),
    "must be a single numeric value"
  )
})

test_that("resolve_log_transform_base errors for base <= 0 or base == 1", {
  expect_error(resolve_log_transform_base(0), "greater than 0")
  expect_error(resolve_log_transform_base(-1), "greater than 0")
  expect_error(resolve_log_transform_base(1), "cannot equal 1")
})

test_that("calc_cpm_df preserves non-integer character rownames", {
  df <- nidap_clean_raw_counts |> as.data.frame()
  rownames(df) <- paste0("row_", seq_len(nrow(df)))
  trans.df <- df
  trans.df[, -1] <- edgeR::cpm(as.matrix(df[, -1]))

  result <- calc_cpm_df(df, feature_id_colname = "Gene")
  expect_equal(rownames(result), rownames(trans.df))
  expect_equal(result, trans.df, ignore_attr = TRUE)
})
