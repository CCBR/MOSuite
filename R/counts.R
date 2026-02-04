#' Calculate counts-per-million (CPM) on raw counts in a multiOmicDataSet
#'
#' @param moo multiOmicDataSet object
#' @param ... additional arguments to pass to edgeR::cpm()
#'
#' @return multiOmicDataSet with cpm-transformed counts
#' @export
#'
#' @examples
#' sample_meta <- data.frame(
#'   sample_id = c("KO_S3", "KO_S4", "WT_S1", "WT_S2"),
#'   condition = factor(
#'     c("knockout", "knockout", "wildtype", "wildtype"),
#'     levels = c("wildtype", "knockout")
#'   )
#' )
#' moo <- create_multiOmicDataSet_from_dataframes(sample_meta, gene_counts) |>
#'   calc_cpm()
#' head(moo@counts$cpm)
calc_cpm <- S7::new_generic("calc_cpm", "moo", function(moo, ...) {
  return(S7::S7_dispatch())
})

S7::method(calc_cpm, multiOmicDataSet) <- function(
  moo,
  feature_id_colname = "gene_id",
  ...
) {
  moo@counts$cpm <- moo@counts$raw |>
    calc_cpm_df(feature_id_colname = feature_id_colname)
  return(moo)
}

#' Calculate CPM on a data frame
#'
#' @inheritParams create_multiOmicDataSet_from_dataframes
#' @param dat data frame of counts with a gene column
#' @param ... additional arguments to pass to edger::cpm()
#'
#' @return cpm-transformed counts as a data frame
#' @keywords internal
#'
calc_cpm_df <- function(dat, feature_id_colname = "gene_id", ...) {
  gene_ids <- dat |> dplyr::pull(feature_id_colname)
  row_names <- rownames(dat)
  dat_cpm <- dat |>
    dplyr::select(-tidyselect::any_of(feature_id_colname)) |>
    as.matrix() |>
    edgeR::cpm(...) |>
    as.data.frame()
  dat_cpm[[feature_id_colname]] <- gene_ids
  rownames(dat_cpm) <- if (
    suppressWarnings(all(!is.na(as.integer(row_names))))
  ) {
    as.integer(row_names)
  } else {
    col
  }
  return(dat_cpm |> dplyr::relocate(tidyselect::all_of(feature_id_colname)))
}

#' Convert a data frame of gene counts to a matrix
#'
#' @inheritParams create_multiOmicDataSet_from_dataframes
#' @param counts_tbl expected feature counts as a dataframe or tibble, with all columns except `feature_id_colname`
#'
#' @return matrix of gene counts with rows as gene IDs
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' counts_dat_to_matrix(head(gene_counts))
#' }
counts_dat_to_matrix <- function(counts_tbl, feature_id_colname = NULL) {
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(counts_tbl)[1]
  }
  counts_dat <- counts_tbl |>
    as.data.frame()
  row.names(counts_dat) <- counts_dat |>
    dplyr::pull(feature_id_colname)
  # convert counts tibble to matrix
  counts_mat <- counts_dat |>
    dplyr::select(-tidyselect::any_of(feature_id_colname)) |>
    as.matrix()
  return(counts_mat)
}

#' Convert all numeric columns in a dataframe to integers
#'
#' Round doubles to integers and convert to integer type
#'
#' @param counts_tbl data frame with numeric columns
#'
#' @return data frame with any numeric columns as integers
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' data.frame(a = c(0, 0.1, 2.3, 5L, 6.9)) |> as_integer_df()
#' }
as_integer_df <- function(counts_tbl) {
  counts_tbl |>
    # deseq2 requires integer counts
    dplyr::mutate(dplyr::across(
      dplyr::where(is.numeric),
      \(x) as.integer(round(x, 0))
    )) |>
    return()
}
