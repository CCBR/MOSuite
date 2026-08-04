write_example_json <- function() {
  j <- list(
    feature_counts_filepath = system.file(
      "extdata",
      "RSEM.genes.expected_count.all_samples.txt.gz",
      package = "MOSuite"
    ),
    sample_meta_filepath = system.file(
      "extdata",
      "sample_metadata.tsv.gz",
      package = "MOSuite"
    ),
    moo_output_rds = "moo.rds"
  )
  return(jsonlite::write_json(j, "tests/testthat/data/example_args.json"))
}

test_that("mosuite cli", {
  command <- paste0(
    system.file("exec", "mosuite", package = "MOSuite"),
    " create_multiOmicDataSet_from_files --json=",
    test_path("data", "example_args.json")
  )
  expect_snapshot(system(command))
})

test_that("cli_exec parses args correctly", {
  expect_equal(cli_exec("do_math"), 3)
  expect_equal(cli_exec(c("do_math", "--subtract", "--no-add")), -1)
  expect_equal(cli_exec(c("do_math", "left=2", "right=3")), 5)
})

test_that("cli_exec --json --debug", {
  expect_equal(
    deparse(cli_exec(
      c(
        "create_multiOmicDataSet_from_files",
        paste0(
          '--json="',
          test_path("data", "example_args.json"),
          '"'
        ),
        "--debug"
      )
    )),
    c(
      paste0(
        "MOSuite::create_multiOmicDataSet_from_files(",
        "feature_counts_filepath = \"inst/extdata/RSEM.genes.expected_count.all_samples.txt.gz\", "
      ),
      "    sample_meta_filepath = \"inst/extdata/sample_metadata.tsv.gz\")"
    )
  )
  expect_error(
    cli_exec(c(
      "filter_counts",
      paste0(
        '--json="',
        test_path("data", "example_args.json"),
        '"'
      ),
      "--debug"
    )),
    "moo_input_rds must be included"
  )
})

test_that("mosuite --help", {
  expect_snapshot(cli_exec("--help"))
  expect_snapshot(system(paste(
    system.file("exec", "mosuite", package = "MOSuite"),
    "--help"
  )))
  expect_snapshot(cli_exec("help"))
  expect_true(inherits(
    cli_exec(c(
      "filter_counts",
      "--help"
    )),
    "help_files_with_topic"
  ))
  expect_error(cli_exec("not_a_function"), "not a known function")
})

test_that("cli_parse handles logical-like strings", {
  expect_true(cli_parse("true"))
  expect_true(cli_parse("True"))
  expect_true(cli_parse("TRUE"))
  expect_false(cli_parse("false"))
  expect_false(cli_parse("False"))
  expect_false(cli_parse("FALSE"))
})

test_that("cli_exec_impl passes positional args", {
  # Positional args are parsed via cli_parse() and appended as unnamed elements.
  expect_equal(cli_exec(c("do_math", "TRUE", "FALSE")), 3)
  expect_equal(cli_exec(c("do_math", "FALSE", "TRUE")), -1)
})

test_that("cli_unknown suggests closest matching function", {
  # Test with a typo that has a close match
  result <- cli_unknown("filter_count", getNamespaceExports("MOSuite"))
  expect_match(result, "filter_count is not a known function")
  expect_match(result, "Did you mean 'filter_counts'")

  # Test with another typo
  result <- cli_unknown("batch_correct_count", getNamespaceExports("MOSuite"))
  expect_match(result, "batch_correct_count is not a known function")
  expect_match(result, "Did you mean 'batch_correct_counts'")

  # Test with completely unrelated name (no suggestions)
  result <- cli_unknown("xyz123", getNamespaceExports("MOSuite"))
  expect_match(result, "xyz123 is not a known function")
  expect_false(grepl("Did you mean", result))
})
