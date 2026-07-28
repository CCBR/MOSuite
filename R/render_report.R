#' Render the template report
#'
#' Copy the Quarto template to the current working directory and render it
#' using `quarto::quarto_render()`. The rendered report will be saved in the
#' current working directory. You can specify additional arguments to
#' `quarto::quarto_render()` to customize the rendering process.
#'
#' You can edit the copy of `report.qmd` in the current working directory to customize the report.
#'
#' @param qmd_template Path to the Quarto report file (default is the template report in the package).
#' @param qmd_src Optional path to copy the Quarto report template to before
#'   rendering. If `NULL` (default), the template will be copied to the current
#'   working directory with the same filename as the template. If a file already
#'   exists at `qmd_src`, it will not be overwritten.
#' @param ... Additional arguments passed to `quarto::quarto_render()`, such as a named list of parameters.
#'
#' @export
#'
#' @examples
#' if (requireNamespace("quarto", quietly = TRUE) &&
#'   requireNamespace("knitr", quietly = TRUE) &&
#'   requireNamespace("rmarkdown", quietly = TRUE)) {
#'   render_report(execute_params = list(
#'     counts_csv = system.file("extdata", "nidap", "Raw_Counts.csv.gz",
#'       package = "MOSuite"
#'     ),
#'     samplesheet_csv = system.file("extdata", "nidap",
#'       "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz",
#'       package = "MOSuite"
#'     )
#'   ))
#' }
render_report <- function(
  qmd_template = system.file("quarto", "report.qmd", package = "MOSuite"),
  qmd_src = NULL,
  ...
) {
  abort_packages_not_installed(c("quarto", "knitr", "rmarkdown")) # nolint: object_usage_linter
  if (is.null(qmd_src)) {
    qmd_src <- basename(qmd_template)
  }
  if (!file.exists(qmd_src)) {
    ok <- file.copy(qmd_template, qmd_src, overwrite = FALSE)
    if (!isTRUE(ok)) {
      stop(glue::glue(
        "Failed to copy template from '{qmd_template}' to '{qmd_src}'"
      ))
    }
  }
  return(quarto::quarto_render(
    input = qmd_src,
    ...
  ))
}
