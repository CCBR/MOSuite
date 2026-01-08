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
  return(jsonlite::write_json(j, "inst/extdata/example.json"))
}

test_that("mosuite cli", {
  command <- paste0(
    system.file("exec", "mosuite", package = "MOSuite"),
    " create_multiOmicDataSet_from_files --json=",
    system.file("extdata", "example.json", package = "MOSuite")
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
          system.file("extdata", "example.json", package = "MOSuite"),
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
        system.file("extdata", "example.json", package = "MOSuite"),
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
  expect_warning(cli_exec("not_a_function"), "not a known function")
})

test_that("mosuite cli E2E", {
  new <- tempfile()
  create_empty_dir(new)
  # note: file paths in json files assume all files are in the current workdir
  withr::with_dir(new = new, code = {
    file.copy(
      system.file("extdata", "nidap", "Raw_Counts.csv.gz", package = "MOSuite"),
      "./"
    )
    file.copy(
      system.file(
        "extdata",
        "nidap",
        "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz",
        package = "MOSuite"
      ),
      "./"
    )
    json_paths <- system.file(
      "extdata",
      "json_args",
      "common",
      package = "MOSuite"
    )
    Sys.glob(glue::glue("{json_paths}/*.json")) |>
      lapply(function(x) {
        file.copy(x, "./")
      })

    run_function_cli("create_multiOmicDataSet_from_files")
    run_function_cli("clean_raw_counts")
    run_function_cli("filter_counts")
    run_function_cli("normalize_counts")
    run_function_cli("batch_correct_counts")
    run_function_cli("diff_counts")
    run_function_cli("filter_diff")

    expect_true(file.exists("moo_diff_filter.rds"))
    moo <- readr::read_rds("moo_diff_filter.rds")
    expect_equal(names(moo@counts), c("raw", "clean", "filt", "norm", "batch"))
    expect_equal(names(moo@analyses), c("colors", "diff", "diff_filt"))
  })
})
