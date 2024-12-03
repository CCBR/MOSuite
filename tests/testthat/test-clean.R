equal_dfs <- function(x, y) {
  all(
    class(x) == class(y),
    names(x) == names(y),
    x == y,
    all.equal(lapply(x, class), lapply(y, class))
  )
}

test_that("clean_raw_counts works", {
  moo <- create_multiOmicDataSet_from_dataframes(
    sample_meta_dat = as.data.frame(nidap_sample_metadata),
    count_dat = as.data.frame(nidap_raw_counts),
    sample_id_colname = "Sample"
  ) %>% clean_raw_counts()
  expect_true(equal_dfs(moo@counts[["clean"]], as.data.frame(nidap_clean_raw_counts)))
})
