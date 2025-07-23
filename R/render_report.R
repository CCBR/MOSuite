#' Render the template report
#'
#' @param input Path to the Quarto report file (default is the template report in the package).
#' @param ... Additional arguments passed to `quarto::quarto_render()`, such as a named list of parameters.
#'
#' @export
#'
#' @example
#' render_report(execute_params = list(counts_tsv = system.file("extdata", "RSEM.genes.expected_count.all_samples.txt.gz", package = "MOSuite"), samplesheet_tsv = system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite")))
render_report <- function(
    input = system.file("quarto", "report.qmd", package = "MOSuite"),
    ...) {
  abort_packages_not_installed(c("quarto", "knitr", "rmarkdown"))
  quarto::quarto_render(
    input = input,
    ...
  )
}
