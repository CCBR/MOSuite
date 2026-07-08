test_that("render_report runs in a temporary directory", {
  skip_if_not_installed("quarto")
  skip_if_not_installed("knitr")
  skip_if_not_installed("rmarkdown")

  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)

  expect_no_error(
    render_report(
      execute_params = list(
        counts_csv = system.file(
          "extdata",
          "nidap",
          "Raw_Counts.csv.gz",
          package = "MOSuite"
        ),
        samplesheet_csv = system.file(
          "extdata",
          "nidap",
          "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz",
          package = "MOSuite"
        ),
        group_colname = "Group",
        label_colname = "Label",
        batch_colname = "Batch",
        contrasts = c("B-A", "C-A", "B-C")
      )
    )
  )

  expect_true(file.exists("report.qmd"))
  expect_true(file.exists("report.html"))
})
