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
#' @param output_dir Directory where the rendered report and supporting files
#'   will be saved. Defaults to `"."` (current working directory). Resolved to
#'   an absolute path before passing to `quarto::quarto_render()`.
#' @param ... Additional arguments passed to `quarto::quarto_render()`, such as a named list of parameters.
#'
#' @export
#'
#' @examples
#' render_report(execute_params = list(
#'   counts_csv = system.file("extdata", "nidap", "Raw_Counts.csv.gz",
#'                            package = "MOSuite"),
#'   samplesheet_csv = system.file("extdata", "nidap",
#'     "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz",
#'     package = "MOSuite")
#' ))
render_report <- function(
  qmd_template = system.file("quarto", "report.qmd", package = "MOSuite"),
  qmd_src = NULL,
  output_dir = ".",
  ...
) {
  abort_packages_not_installed(c("quarto", "knitr", "rmarkdown")) # nolint: object_usage_linter
  if (is.null(qmd_src)) {
    qmd_src <- basename(qmd_template)
  }
  # Resolve to absolute path so quarto subprocess uses the correct directory
  qmd_src <- normalizePath(qmd_src, mustWork = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
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
    output_dir = output_dir,
    ...
  ))
}
