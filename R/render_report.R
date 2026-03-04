#' Render the template report
#'
#' Copy the Quarto template to the current working directory and render it using `quarto::quarto_render()`. The rendered report will be saved in the current working directory. You can specify additional arguments to `quarto::quarto_render()` to customize the rendering process.
#'
#' @param qmd_template Path to the Quarto report file (default is the template report in the package).
#' @param ... Additional arguments passed to `quarto::quarto_render()`, such as a named list of parameters.
#'
#' @export
#'
#' @example
#' render_report(execute_params = list(counts_tsv = system.file("extdata", "RSEM.genes.expected_count.all_samples.txt.gz", package = "MOSuite"), samplesheet_tsv = system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite")))
render_report <- function(
  qmd_template = system.file("quarto", "report.qmd", package = "MOSuite"),
  ...
) {
  abort_packages_not_installed(c("quarto", "knitr", "rmarkdown"))
  qmd_src = basename(qmd_template)
  file.copy(qmd_template, qmd_src, overwrite = TRUE)
  return(quarto::quarto_render(
    input = qmd_src,
    ...
  ))
}
