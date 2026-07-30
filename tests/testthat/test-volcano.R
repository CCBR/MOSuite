df_both <- data.frame(
  Gene = letters[1:3],
  `B-A_logFC` = c(1, -1, 0.1),
  `B-A_adjpval` = c(0.01, 0.02, 0.5),
  `B-A_pval` = c(0.001, 0.002, 0.4),
  `C-A_logFC` = c(2, -2, 0.2),
  `C-A_adjpval` = c(0.03, 0.04, 0.6),
  `C-A_pval` = c(0.003, 0.004, 0.5),
  check.names = FALSE
)

df_pval_only <- data.frame(
  Gene = letters[1:3],
  `B-A_logFC` = c(1, -1, 0.1),
  `B-A_pval` = c(0.001, 0.002, 0.4),
  check.names = FALSE
)

df_no_signif <- data.frame(
  Gene = letters[1:3],
  `B-A_logFC` = c(1, -1, 0.1),
  check.names = FALSE
)

df_no_logfc <- data.frame(
  Gene = letters[1:3],
  `B-A_adjpval` = c(0.01, 0.02, 0.5),
  check.names = FALSE
)

test_that("resolve_volcano_colnames auto-detects both when NULL", {
  result <- resolve_volcano_colnames(df_both, NULL, NULL)
  expect_equal(result$change_colname, c("B-A_logFC", "C-A_logFC"))
  expect_equal(result$signif_colname, c("B-A_adjpval", "C-A_adjpval"))
})

test_that("resolve_volcano_colnames prefers adjpval over pval", {
  result <- resolve_volcano_colnames(df_both, NULL, NULL)
  expect_true(all(grepl("_adjpval$", result$signif_colname)))
})

test_that("resolve_volcano_colnames falls back to pval when no adjpval", {
  result <- resolve_volcano_colnames(df_pval_only, NULL, NULL)
  expect_equal(result$change_colname, "B-A_logFC")
  expect_equal(result$signif_colname, "B-A_pval")
})

test_that("resolve_volcano_colnames respects explicit change_colname", {
  result <- resolve_volcano_colnames(df_both, "C-A_logFC", NULL)
  expect_equal(result$change_colname, "C-A_logFC")
  expect_equal(result$signif_colname, "C-A_adjpval")
})

test_that("resolve_volcano_colnames passes through explicit signif_colname", {
  result <- resolve_volcano_colnames(
    df_both,
    c("B-A_logFC", "C-A_logFC"),
    c("B-A_pval", "C-A_pval")
  )
  expect_equal(result$change_colname, c("B-A_logFC", "C-A_logFC"))
  expect_equal(result$signif_colname, c("B-A_pval", "C-A_pval"))
})

test_that("resolve_volcano_colnames errors when no _logFC columns exist", {
  expect_error(
    resolve_volcano_colnames(df_no_logfc, NULL, NULL),
    regexp = "_logFC"
  )
})

test_that("resolve_volcano_colnames errors when no significance columns found", {
  expect_error(
    resolve_volcano_colnames(df_no_signif, NULL, NULL),
    regexp = "auto-detect"
  )
})

test_that("resolve_volcano_colnames errors when lengths differ", {
  expect_error(
    resolve_volcano_colnames(
      df_both,
      c("B-A_logFC", "C-A_logFC"),
      "B-A_adjpval"
    ),
    regexp = "same length"
  )
})
