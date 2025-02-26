#' Plot read depth as a bar plot
#'
#' The first argument can be a `multiOmicDataset` object (`moo`) or a `data.frame` containing counts.
#' For a `moo`, choose which counts slot to use with `count_type` & `sub_count_type`.
#'
#' # Methods
#'
#'  See documentation below for method-specific arguments
#'
#'   | method                   | class of `moo_counts` |
#'   |--------------------------|--------------------|
#'   | [plot_read_depth_moo]     | `multiOmicDataset` |
#'   | [plot_read_depth_dat]     | `data.frame`       |
#'
#' @export
#' @examples
#' moo <- multiOmicDataSet(
#'   sample_metadata = nidap_sample_metadata,
#'   anno_dat = data.frame(),
#'   counts_lst = list("raw" = nidap_raw_counts)
#' )
#'
#' # multiOmicDataSet
#' plot_read_depth(moo)
#'
#' # dataframe
#' plot_read_depth(nidap_clean_raw_counts)
#'
#' @family plotters
plot_read_depth <- S7::new_generic("plot_read_depth", dispatch_args = "moo_counts")

#' Plot read depth for multiOmicDataSet
#'
#' @param moo_counts multiOmicDataSet containing `count_type` & `sub_count_type` in the counts slot
#' @param count_type the type of counts to use -- must be a name in the counts slot (`moo@counts[[count_type]]`)
#' @param sub_count_type if `count_type` is a list, specify the sub count type within the list (`moo@counts[[count_type]][[sub_count_type]]`)
#' @param ... remaining arguments forwarded to the `plot_read_depth` method for `data.frame`s (see `?plot_read_depth_dat`)
#'
#' @returns ggplot object
#' @examples
#' moo <- multiOmicDataSet(
#'   sample_metadata = nidap_sample_metadata,
#'   anno_dat = data.frame(),
#'   counts_lst = list("raw" = nidap_raw_counts)
#' )
#'
#' plot_read_depth(moo)
#' @name plot_read_depth_moo
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
#' @param moo_counts data.frame with counts data, e.g. from the counts slot of a [multiOmicDataSet]
#'
#' @returns ggplot object
#' @examples
#' plot_read_depth(nidap_clean_raw_counts)
#'
#' @name plot_read_depth_dat
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
