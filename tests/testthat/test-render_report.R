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
        ),
        group_colname = 'condition',
        label_colname = 'sample_id',
        batch_colname = 'sample_id',
        contrasts = c("knockout-wildtype")
      )
    )
  )

  expect_true(file.exists("report.qmd"))
  expect_true(file.exists("report.html"))

  expect_no_error(
    render_report(
      execute_params = list(
        counts_tsv = system.file(
          "extdata",
          "nidap",
          "Raw_Counts.csv.gz",
          package = "MOSuite"
        ),
        samplesheet_tsv = system.file(
          "extdata",
          "nidap",
          "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz",
          package = "MOSuite"
        )
      )
    )
  )
})
