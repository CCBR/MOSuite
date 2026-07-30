#' Resolve volcano plot column names
#'
#' Auto-detects `change_colname` and `signif_colname` from a data frame when
#' either is `NULL`. Used by [plot_volcano_enhanced()] and
#' [plot_volcano_summary()].
#'
#' When `change_colname` is `NULL`, all columns ending in `_logFC` are used.
#' When `signif_colname` is `NULL`, significance columns are detected by
#' checking for `_adjpval` columns first, then `_pval`, for each contrast
#' derived from `change_colname`.
#'
#' @param diff_dat A data frame of differential analysis results.
#' @param change_colname Character vector of logFC column names, or `NULL` to
#'   auto-detect columns ending in `_logFC`.
#' @param signif_colname Character vector of significance column names, or
#'   `NULL` to auto-detect, preferring `_adjpval` over `_pval`.
#' @return A named list with elements `change_colname` and `signif_colname`.
#' @keywords internal
resolve_volcano_colnames <- function(diff_dat, change_colname, signif_colname) {
  if (is.null(change_colname)) {
    change_colname <- grep("_logFC$", colnames(diff_dat), value = TRUE)
    if (length(change_colname) == 0) {
      cli::cli_abort(
        "No columns ending in {.val _logFC} found. \\
        Supply {.arg change_colname} explicitly."
      )
    }
  }
  if (is.null(signif_colname)) {
    contrasts <- sub("_logFC$", "", change_colname)
    detected_suffix <- NULL
    for (suffix in c("adjpval", "pval")) {
      if (all(paste0(contrasts, "_", suffix) %in% colnames(diff_dat))) {
        detected_suffix <- suffix
        break
      }
    }
    if (is.null(detected_suffix)) {
      cli::cli_abort(
        "Could not auto-detect significance columns. \\
        Supply {.arg signif_colname} explicitly \\
        (e.g. {.code c(\"B-A_adjpval\", \"C-A_adjpval\")})."
      )
    }
    signif_colname <- paste0(contrasts, "_", detected_suffix)
  }
  if (length(change_colname) != length(signif_colname)) {
    cli::cli_abort(
      "{.arg change_colname} (length {length(change_colname)}) and \\
      {.arg signif_colname} (length {length(signif_colname)}) must have the \\
      same length."
    )
  }
  list(change_colname = change_colname, signif_colname = signif_colname)
}
