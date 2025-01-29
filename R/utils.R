#' Glue gene_id and GeneName columns into one column
#'
#' @param count_dat data frame containing gene_id and GeneName columns
#'
#' @returns count_dat with gene_id and GeneName joined with `|` as the new gene_id column
#' @keywords internal
#' @examples
#' \dontrun{
#' gene_counts %>%
#'   glue_gene_symbols() %>%
#'   head()
#' }
glue_gene_symbols <- function(count_dat) {
  if ("gene_id" %in% colnames(count_dat) & "GeneName" %in% colnames(count_dat)) {
    count_dat <- count_dat %>%
      dplyr::mutate(
        gene_id = glue::glue("{gene_id}|{GeneName}"),
        .keep = "unused"
      )
  }
  return(count_dat)
}
