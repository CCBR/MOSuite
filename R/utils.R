#' Glue gene_id and GeneName columns into one column
#'
#' @param counts_dat data frame containing gene_id and GeneName columns
#'
#' @returns counts_dat with gene_id and GeneName joined with `|` as the new gene_id column
#' @keywords internal
#' @examples
#' \dontrun{
#' gene_counts %>%
#'   glue_gene_symbols() %>%
#'   head()
#' }
glue_gene_symbols <- function(counts_dat) {
  if ("gene_id" %in% colnames(counts_dat) & "GeneName" %in% colnames(counts_dat)) {
    counts_dat <- counts_dat %>%
      dplyr::mutate(
        gene_id = glue::glue("{gene_id}|{GeneName}"),
        .keep = "unused"
      )
  }
  return(counts_dat)
}

#' Check whether package(s) are installed
#'
#' @param ... names of packages to check
#' @return named vector with status of each packages; installed (`TRUE`) or not (`FALSE`)
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' check_packages_installed("base")
#' check_packages_installed("not-a-package-name")
#' all(check_packages_installed("parallel", "doFuture"))
#' }
check_packages_installed <- function(...) {
  return(sapply(c(...), requireNamespace, quietly = TRUE))
}

#' Throw error if required packages are not installed.
#'
#' Reports which packages need to be installed and the parent function name.
#' See \url{https://stackoverflow.com/questions/15595478/how-to-get-the-name-of-the-calling-function-inside-the-called-routine}
#' This is only intended to be used inside a function. It will error otherwise.
#'
#' @inheritParams check_packages_installed
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' abort_packages_not_installed("base")
#' abort_packages_not_installed("not-a-package-name", "caret", "dplyr", "non_package")
#' }
abort_packages_not_installed <- function(...) {
  package_status <- check_packages_installed(...)
  parent_fcn_name <- sub("\\(.*$", "\\(\\)", deparse(sys.calls()[[sys.nframe() - 1]]))
  packages_not_installed <- Filter(isFALSE, package_status)
  if (length(packages_not_installed) > 0) {
    msg <- paste0(
      "The following package(s) are required for `", parent_fcn_name,
      "` but are not installed: \n  ",
      paste0(names(packages_not_installed), collapse = ", ")
    )
    stop(msg)
  }
}
