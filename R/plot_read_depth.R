#' Plot read depth as a bar plot
#'
#' The first argument can be a `multiOmicDataset` object (`moo`) or a `data.frame` containing counts.
#' For a `moo`, choose which counts slot to use with `count_type` & (optionally) `sub_count_type`.
#'
#' @param moo_counts counts dataframe or `multiOmicDataSet` containing `count_type` & `sub_count_type` in the counts slot
#' @param ... arguments forwarded to method
#'
#' @return ggplot barplot
#'
#' @export
#' @examples
#' # multiOmicDataSet
#' moo <- multiOmicDataSet(
#'   sample_metadata = nidap_sample_metadata,
#'   anno_dat = data.frame(),
#'   counts_lst = list(
#'     "raw" = nidap_raw_counts,
#'     "clean" = nidap_clean_raw_counts
#'   )
#' )
#'
#' plot_read_depth(moo, count_type = "clean")
#'
#' # dataframe
#' plot_read_depth(nidap_clean_raw_counts)
#'
#' @details
#'
#' # Methods
#'
#' | link to docs  | class  |
#' |---|---|
#' | [plot_read_depth_moo] | `multiOmicDataSet` |
#' | [plot_read_depth_dat] | `data.frame`       |
#'
#' @family plotters
#' @keywords plotters
#' @family moo methods
plot_read_depth <- S7::new_generic("plot_read_depth", dispatch_args = "moo_counts")

#' Plot read depth for multiOmicDataSet
#'
#' @name plot_read_depth_moo
#'
#' @param moo_counts `multiOmicDataSet` containing `count_type` & `sub_count_type` in the counts slot
#' @param count_type the type of counts to use. Must be a name in the counts slot (`names(moo@counts)`).
#' @param sub_count_type used if `count_type` is a list in the counts slot: specify the sub count type within the list. Must be a name in `names(moo@counts[[count_type]])`.
#'
#' @return ggplot barplot
#'
#' @examples
#' # multiOmicDataSet
#' moo <- multiOmicDataSet(
#'   sample_metadata = nidap_sample_metadata,
#'   anno_dat = data.frame(),
#'   counts_lst = list(
#'     "raw" = nidap_raw_counts,
#'     "clean" = nidap_clean_raw_counts
#'   )
#' )
#'
#' plot_read_depth(moo, count_type = "clean")
#'
#' @seealso [plot_read_depth] generic
#' @family plotters for multiOmicDataSets
S7::method(plot_read_depth, multiOmicDataSet) <- function(moo_counts,
                                                          count_type,
                                                          sub_count_type = NULL,
                                                          ...) {
  counts_dat <- extract_counts(moo_counts, count_type, sub_count_type)
  plot_read_depth(counts_dat, ...)
}

#' Plot read depth for `data.frame`
#'
#' @name plot_read_depth_dat
#'
#' @param moo_counts counts dataframe
#'
#' @return ggplot barplot
#'
#' @examples
#' # dataframe
#' plot_read_depth(nidap_clean_raw_counts)
#'
#' @seealso [plot_read_depth] generic
#' @family plotters for counts dataframes
S7::method(plot_read_depth, S7::class_data.frame) <- function(moo_counts) {
  sample_names <- read_sums <- column_sums <- NULL
  counts_dat <- moo_counts
  sum_df <- counts_dat %>%
    dplyr::summarize(dplyr::across(tidyselect::where(is.numeric), sum)) %>%
    tidyr::pivot_longer(dplyr::everything(),
      names_to = "sample_names",
      values_to = "column_sums"
    )

  # Plotting
  read_plot <- ggplot2::ggplot(sum_df, ggplot2::aes(x = sample_names, y = column_sums)) +
    ggplot2::geom_bar(stat = "identity", fill = "blue") +
    ggplot2::labs(title = "Total Reads per Sample", x = "Samples", y = "Read Count") +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 14
      ),
      axis.text.y = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 16),
      plot.title = ggplot2::element_text(size = 20)
    )
  return(read_plot)
}
