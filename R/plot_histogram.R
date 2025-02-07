#' Plot histogram
#'
#' @inheritParams filter_counts
#' @param counts_dat counts data
#' @param sample_metadata sample metadata as a data frame or tibble.
#' @param x_axis_label text label for the x axis `ggplot2::xlab()`
#' @param y_axis_label text label for the y axis `ggplot2::ylab()`
#' @param color_values vector of colors as hex values or names recognized by R
#' @param color_by_group Set to FALSE to label histogram by Sample Names, or set to TRUE to label histogram by the column you select in the "Group Column Used to Color Histogram" parameter (below). Default is FALSE.
#' @param set_min_max_for_x_axis whether to override the default for `ggplot2::xlim()` (default: `FALSE`)
#' @param minimum_for_x_axis value to override default `min` for `ggplot2::xlim()`
#' @param maximum_for_x_axis value to override default `max` for `ggplot2::xlim()`
#' @param legend_position passed to in `legend.position` `ggplot2::theme()`
#' @param legend_font_size passed to `ggplot2::element_text()` via `ggplot2::theme()`
#' @param number_of_legend_columns passed to `ncol` in `ggplot2::guide_legend()`
#'
#' @returns ggplot object
#' @export
#'
plot_histogram <- function(counts_dat,
                           sample_metadata,
                           sample_id_colname = NULL,
                           feature_id_colname = NULL,
                           group_colname = "Group",
                           label_colname = "Label",
                           color_values = c(
                             "#5954d6", "#e1562c", "#b80058", "#00c6f8", "#d163e6", "#00a76c",
                             "#ff9287", "#008cf9", "#006e00", "#796880", "#FFA500", "#878500"
                           ),
                           color_by_group = FALSE,
                           set_min_max_for_x_axis = FALSE,
                           minimum_for_x_axis = -1,
                           maximum_for_x_axis = 1,
                           legend_position = "top",
                           legend_font_size = 10,
                           number_of_legend_columns = 6,
                           x_axis_label = "Counts",
                           y_axis_label = "Density",
                           interactive_plots = FALSE) {
  count <- NULL

  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(counts_dat)[1]
  }

  df_long <- counts_dat %>%
    tidyr::pivot_longer(-tidyselect::all_of(feature_id_colname),
      names_to = sample_id_colname,
      values_to = "count"
    ) %>%
    dplyr::left_join(sample_metadata, by = sample_id_colname)

  if (set_min_max_for_x_axis == TRUE) {
    xmin <- minimum_for_x_axis
    xmax <- maximum_for_x_axis
  } else {
    xmin <- min(df_long %>% dplyr::pull(count))
    xmax <- max(df_long %>% dplyr::pull(count))
  }

  if (color_by_group == TRUE) {
    df_long %<>%
      dplyr::mutate(!!rlang::sym(group_colname) := as.factor(!!rlang::sym(group_colname))) %>%
      dplyr::filter(!is.na(group_colname))
    n <- df_long %>%
      dplyr::pull(group_colname) %>%
      levels() %>%
      length()

    # plot Density
    hist_plot <- df_long %>%
      ggplot2::ggplot(ggplot2::aes(x = count, group = !!rlang::sym(sample_id_colname))) +
      ggplot2::geom_density(ggplot2::aes(colour = !!rlang::sym(group_colname)), linewidth = 1)
  } else {
    n <- df_long %>%
      dplyr::pull(sample_id_colname) %>%
      unique() %>%
      length()

    hist_plot <- df_long %>%
      ggplot2::ggplot(ggplot2::aes(x = count, group = !!rlang::sym(sample_id_colname))) +
      ggplot2::geom_density(ggplot2::aes(colour = !!rlang::sym(sample_id_colname)), linewidth = 1)
  }

  hist_plot <- hist_plot +
    ggplot2::xlab(x_axis_label) +
    ggplot2::ylab(y_axis_label) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = legend_position,
      legend.text = ggplot2::element_text(size = legend_font_size),
      legend.title = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 18),
      axis.title = ggplot2::element_text(size = 20),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0),
      axis.line = ggplot2::element_line(linewidth = .5),
      axis.ticks = ggplot2::element_line(linewidth = 1)
    ) +
    ggplot2::ggtitle("Frequency Histogram") +
    ggplot2::xlim(xmin, xmax) +
    # scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),n)) +
    ggplot2::scale_colour_manual(values = color_values[1:n]) +
    ggplot2::guides(linetype = ggplot2::guide_legend(ncol = number_of_legend_columns))

  if (isTRUE(interactive_plots)) {
    hist_plot <- (hist_plot + ggplot2::theme(legend.position = "none")) %>%
      plotly::ggplotly(tooltip = c(sample_id_colname))
  }
  return(hist_plot)
}
