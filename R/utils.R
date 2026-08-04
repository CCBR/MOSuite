#' Glue gene_id and GeneName columns into one column
#'
#' @param counts_dat data frame containing gene_id and GeneName columns
#'
#' @returns counts_dat with gene_id and GeneName joined with `|` as the new gene_id column
#' @keywords internal
#' @examples
#' \dontrun{
#' gene_counts |>
#'   glue_gene_symbols() |>
#'   head()
#' }
glue_gene_symbols <- function(counts_dat) {
  if (
    "gene_id" %in% colnames(counts_dat) && "GeneName" %in% colnames(counts_dat)
  ) {
    counts_dat <- counts_dat |>
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
  packages_not_installed <- Filter(isFALSE, package_status)
  if (length(packages_not_installed) > 0) {
    msg <- paste0(
      "The following package(s) are required but are not installed: \n  ",
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
do_math <- function(add = TRUE, subtract = FALSE, left = 1, right = 2) {
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
#' dfs |> join_dfs_wide()
#'
join_dfs_wide <- function(df_list, join_fn = dplyr::left_join) {
  if (!inherits(df_list, "list")) {
    stop(glue::glue("df_list must be a named list. class: {class(df_list)}"))
  }
  if (is.null(names(df_list))) {
    stop(glue::glue("df_list does not have names"))
  }
  # use first column as start
  common_col <- df_list[[1]] |>
    dplyr::select(1) |>
    colnames()
  dat_joined <- purrr::map(names(df_list), \(df_name) {
    result <- df_list[[df_name]] |>
      dplyr::rename_with(
        .cols = !tidyselect::any_of(common_col),
        .fn = \(x) {
          return(glue::glue("{df_name}_{x}"))
        }
      )
    return(result)
  }) |>
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
#' dfs |> bind_dfs_long()
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
  common_col <- df_list[[1]] |>
    dplyr::select(1) |>
    colnames()
  dat_joined <- purrr::map(names(df_list), \(df_name) {
    result <- df_list[[df_name]] |>
      dplyr::mutate({{ outcolname }} := df_name, .after = common_col)
    return(result)
  }) |>
    dplyr::bind_rows()
  return(dat_joined)
}

#' Set up capsule environment and directories
#'
#' Initializes the results directory structure and logs installed R package versions.
#' This is a common setup task used across all Code Ocean capsules.
#'
#' @param base_results_dir base path to results directory (default: `../results`)
#'
#' @returns invisibly returns a list with `results_dir` and `plots_dir` paths
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' setup_capsule_environment()
#' }
#'
#' @export
setup_capsule_environment <- function(
  base_results_dir = file.path("..", "results")
) {
  results_dir <- base_results_dir
  plots_dir <- file.path(results_dir, "figures")

  # Set options for plots directory
  options(moo_plots_dir = plots_dir, moo_save_plots = TRUE)

  # Log installed packages & versions
  pkg_versions <- tibble::as_tibble(utils::installed.packages())
  readr::write_csv(pkg_versions, file.path(results_dir, "r-packages.csv"))

  return(invisible(list(results_dir = results_dir, plots_dir = plots_dir)))
}

#' Load multiOmicDataSet from data directory
#'
#' Searches the ../data directory for .rds files and loads the first matching
#' multiOmicDataSet object. Validates that the loaded object is of the correct class.
#'
#' @param data_dir path to data directory containing .rds file (default: `../data`)
#'
#' @returns loaded multiOmicDataSet object
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' moo <- load_moo_from_data_dir()
#' }
#'
#' @export
load_moo_from_data_dir <- function(data_dir = file.path("..", "data")) {
  regex_moo <- ".*\\.rds$"
  data_files <- list.files(data_dir, recursive = TRUE, full.names = TRUE)
  moo_files <- Filter(
    \(x) stringr::str_detect(x, stringr::regex(regex_moo, ignore_case = TRUE)),
    data_files
  )

  if (length(moo_files) == 0) {
    stop(glue::glue("No files matching regex: {regex_moo}"))
  }

  moo_filename <- moo_files[1]
  moo <- readr::read_rds(moo_filename)

  message(glue::glue("Reading multiOmicDataSet from {moo_filename}"))

  if (!inherits(moo, "MOSuite::multiOmicDataSet")) {
    stop(glue::glue("The input is not a multiOmicDataSet. class: {class(moo)}"))
  }

  return(moo)
}

#' Parse comma-separated string into a vector
#'
#' Splits a comma-separated string into a trimmed character vector.
#' Returns NULL if input is empty, NULL, or has zero length.
#'
#' @param x character string with comma-separated values
#'
#' @returns character vector or NULL if input is empty
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' parse_optional_vector("a, b, c")
#' parse_optional_vector("")
#' }
#'
#' @export
parse_optional_vector <- function(x) {
  if (is.null(x) || identical(x, "") || length(x) == 0) {
    return(NULL)
  }
  return(trimws(unlist(strsplit(x, ","))))
}

#' Parse comma-separated string with default fallback
#'
#' Splits a comma-separated string into a trimmed character vector.
#' Returns a default value if input is empty, NULL, or has zero length.
#'
#' @param x character string with comma-separated values
#' @param default default value to return if x is empty
#'
#' @returns character vector or default value
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' parse_vector_with_default("a, b, c", "default")
#' parse_vector_with_default("", "default")
#' }
#'
#' @export
parse_vector_with_default <- function(x, default) {
  parsed <- parse_optional_vector(x)
  if (is.null(parsed)) {
    return(default)
  }
  return(parsed)
}

#' Parse sample rename pairs from string
#'
#' Parses a string containing sample rename pairs in format "old:new,old2:new2"
#' and returns a named list where names are old sample names and values are new names.
#'
#' @param x character string with rename pairs in format "old:new,old2:new2"
#'
#' @returns named list with old names as keys and new names as values, or NULL if empty
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' parse_samples_to_rename("sample1:S1,sample2:S2")
#' parse_samples_to_rename("")
#' }
#'
#' @export
parse_samples_to_rename <- function(x) {
  if (is.null(x) || identical(x, "") || length(x) == 0) {
    return(NULL)
  }

  pairs <- trimws(unlist(strsplit(x, ",")))
  result <- list()

  for (pair in pairs) {
    parts <- trimws(unlist(strsplit(pair, ":")))
    if (length(parts) == 2) {
      result[[parts[1]]] <- parts[2]
    }
  }

  if (length(result) == 0) {
    return(NULL)
  }

  return(result)
}
