test_that("render_report runs in a temporary directory", {
  skip_if_not_installed("quarto")

  withr::local_dir(tempdir())

  expect_no_error(
    render_report(
      execute_params = list(
        counts_tsv = system.file(
          "extdata",
          "RSEM.genes.expected_count.all_samples.txt.gz",
          package = "MOSuite"
        ),
        samplesheet_tsv = system.file(
          "extdata",
          "sample_metadata.tsv.gz",
          package = "MOSuite"
        )
      )
    )
  )

  expect_true(file.exists("report.qmd"))
  expect_true(file.exists("report.html"))
})
