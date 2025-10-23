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
  if ("gene_id" %in% colnames(counts_dat) && "GeneName" %in% colnames(counts_dat)) {
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
#' See
#' https://stackoverflow.com/questions/15595478/how-to-get-the-name-of-the-calling-function-inside-the-called-routine
#'
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
      "The following package(s) are required for `",
      parent_fcn_name,
      "` but are not installed: \n  ",
      paste0(names(packages_not_installed), collapse = ", ")
    )
    stop(msg)
  }
}

#' Function for testing CLI argument parsing
#'
#' @param add whether to add left and right
#' @param subtract whether to subtract left and right
#' @param left number on the left side of the operand
#' @param right number on the right side of the operand
#'
#' @returns result of adding or subtracting left and right
#'
#' @export
#' @keywords internal
do_math <- function(add = TRUE,
                    subtract = FALSE,
                    left = 1,
                    right = 2) {
  result <- NULL
  if (isTRUE(add)) {
    result <- left + right
  } else if (isTRUE(subtract)) {
    result <- left - right
  }
  return(result)
}

#' Join dataframes in named list to wide dataframe
#'
#' The first column is assumed to be shared by all dataframes
#'
#' @param df_list named list of dataframes
#' @param join_fn join function to use (Default: `dplyr::left_join`)
#'
#' @returns wide dataframe
#' @export
#' @keywords utilities
#'
#' @examples
#'
#' dfs <- list(
#'   "a_vs_b" = data.frame(id = c("a1", "b2", "c3"), score = runif(3)),
#'   "b_vs_c" = data.frame(id = c("a1", "b2", "c3"), score = rnorm(3))
#' )
#' dfs %>% join_dfs_wide()
#'
join_dfs_wide <- function(df_list, join_fn = dplyr::left_join) {
  if (!inherits(df_list, "list")) {
    stop(glue::glue("df_list must be a named list. class: {class(df_list)}"))
  }
  if (is.null(names(df_list))) {
    stop(glue::glue("df_list does not have names"))
  }
  # use first column as start
  common_col <- df_list[[1]] %>%
    dplyr::select(1) %>%
    colnames()
  dat_joined <- purrr::map(names(df_list), \(df_name) {
    df_list[[df_name]] %>%
      dplyr::rename_with(
        .cols = !tidyselect::any_of(common_col),
        .fn = \(x) {
          return(glue::glue("{df_name}_{x}"))
        }
      ) %>%
      return()
  }) %>%
    purrr::reduce(join_fn)
  return(dat_joined)
}

#' Bind dataframes in named list to long dataframe
#'
#' The dataframes must have all of the same columns
#'
#' @param df_list named list of dataframes
#' @param outcolname column name in output dataframe for the names from the named list
#'
#' @returns long dataframe with new column `outcolname` from named list
#' @export
#' @keywords utilities
#'
#' @examples
#'
#' dfs <- list(
#'   "a_vs_b" = data.frame(id = c("a1", "b2", "c3"), score = runif(3)),
#'   "b_vs_c" = data.frame(id = c("a1", "b2", "c3"), score = rnorm(3))
#' )
#' dfs %>% bind_dfs_long()
#'
bind_dfs_long <- function(df_list, outcolname = contrast) {
  contrast <- NULL # data variable
  if (!inherits(df_list, "list")) {
    stop(glue::glue("df_list must be a named list. class: {class(df_list)}"))
  }
  if (is.null(names(df_list))) {
    stop(glue::glue("df_list does not have names"))
  }
  # use first column as start
  common_col <- df_list[[1]] %>%
    dplyr::select(1) %>%
    colnames()
  dat_joined <- purrr::map(names(df_list), \(df_name) {
    df_list[[df_name]] %>%
      dplyr::mutate({{ outcolname }} := df_name, .after = common_col) %>%
      return()
  }) %>%
    dplyr::bind_rows()
  return(dat_joined)
}
