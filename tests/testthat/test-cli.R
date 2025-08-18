write_example_json <- function() {
  j <- list(
    feature_counts_filepath = system.file("extdata", "RSEM.genes.expected_count.all_samples.txt.gz", package = "MOSuite"),
    sample_meta_filepath = system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite"),
    moo_output_rds = "moo.rds"
  )
  jsonlite::write_json(j, "inst/extdata/example.json")
}

test_that("mosuite cli", {
  command <- paste0(
    system.file("exec", "mosuite", package = "MOSuite"),
    " create_multiOmicDataSet_from_files --json=",
    system.file("extdata", "example.json", package = "MOSuite")
  )
  expect_snapshot(system(command))
})

test_that("cli_exec --json --debug", {
  expect_snapshot(
    print(cli_exec(c(
      "create_multiOmicDataSet_from_files",
      paste0('--json="', system.file("extdata", "example.json", package = "MOSuite"), '"'),
      "--debug"
    )))
  )
})

test_that("mosuite --help", {
  expect_snapshot(cli_exec("--help"))
  expect_snapshot(system(paste(system.file("exec", "mosuite", package = "MOSuite"), "--help")))
  expect_snapshot(cli_exec(c("filter_counts", "--help")))
  expect_warning(cli_exec("not_a_function"), "not a known function")
})
