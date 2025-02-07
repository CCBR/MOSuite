#' Create read depth plot
#'
#' @param counts_dat dataframe with raw counts data
#'
#' @returns ggplot object
#' @export
#'
plot_read_depth <- function(counts_dat) {
  sample_names <- read_sums <- NULL

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
