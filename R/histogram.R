#' Plot histogram
#'
#' @inheritParams filter_counts
#' @param log_counts log-transformed filtered counts
#' @param sample_metadata sample metadata as a data frame or tibble.
#' @param x_axis_label text label for the x axis
#' @param y_axis_label text label for the y axis
#' @param color_values vector of colors as hex values or names recognized by R
#'
#' @returns ggplot object
#' @export
#'
plot_histogram <- function(log_counts,
                           sample_metadata,
                           feature_id_colname,
                           group_colname,
                           label_colname,
                           color_values,
                           color_histogram_by_group = FALSE,
                           set_min_max_for_x_axis_for_histogram = FALSE,
                           minimum_for_x_axis_for_histogram = -1,
                           maximum_for_x_axis_for_histogram = 1,
                           legend_position_for_histogram = "top",
                           legend_font_size_for_histogram = 10,
                           number_of_histogram_legend_columns = 6,
                           x_axis_label = "Counts",
                           y_axis_label = "Density",
                           make_plots_interactive = FALSE) {
  Var2 <- colgroup <- value <- NULL
  df.m <- reshape2::melt(log_counts, id.vars = c(feature_id_colname))
  df.m <- dplyr::rename(df.m, sample = Var2)

  if (set_min_max_for_x_axis_for_histogram == TRUE) {
    xmin <- minimum_for_x_axis_for_histogram
    xmax <- maximum_for_x_axis_for_histogram
  } else {
    xmin <- min(df.m$value)
    xmax <- max(df.m$value)
  }

  if (color_histogram_by_group == TRUE) {
    df.m <- df.m %>% dplyr::mutate(colgroup = sample_metadata[sample, group_colname])
    df.m <- df.m[stats::complete.cases(df.m[, "colgroup"]), ]
    df.m$colgroup <- gsub("\\s", "_", df.m$colgroup)
    df.m$colgroup <- factor(df.m$colgroup, levels = unique(df.m$colgroup))
    ## print(unique(df.m$sample))
    n <- length(levels(df.m$colgroup))
    cols <- color_values[1:n]

    # plot Density
    hist_plot <- df.m %>%
      ggplot2::ggplot(ggplot2::aes(x = value, group = sample)) +
      ggplot2::geom_density(ggplot2::aes(colour = colgroup), linewidth = 1)
  } else {
    df.m$sample <- sample_metadata[df.m$sample, label_colname]
    n <- length(unique(df.m$sample))
    cols <- get_random_colors(n)

    hist_plot <- df.m %>%
      ggplot2::ggplot(ggplot2::aes(x = value, group = sample)) +
      ggplot2::geom_density(ggplot2::aes(colour = sample), linewidth = 1)
  }

  hist_plot <- hist_plot +
    ggplot2::xlab(x_axis_label) +
    ggplot2::ylab(y_axis_label) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = legend_position_for_histogram,
      legend.text = ggplot2::element_text(size = legend_font_size_for_histogram),
      legend.title = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 18),
      axis.title = ggplot2::element_text(size = 20),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0),
      axis.line = ggplot2::element_line(linewidth = .5),
      axis.ticks = ggplot2::element_line(linewidth = 1)
    ) +
    # ggtitle("Frequency Histogram") +
    ggplot2::xlim(xmin, xmax) +
    # scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),n)) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::guides(linetype = ggplot2::guide_legend(ncol = number_of_histogram_legend_columns))

  if (isTRUE(make_plots_interactive)) {
    hist_plot <- (hist_plot + ggplot2::theme(legend.position = "none")) %>% plotly::ggplotly(tooltip = c("sample"))
  }
  return(hist_plot)
}
