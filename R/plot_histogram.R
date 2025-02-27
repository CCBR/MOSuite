#' Plot histogram
#'
#' @inherit moo_counts description
#'
#'
#' @param moo_counts counts dataframe or `multiOmicDataSet` containing `count_type` & `sub_count_type` in the counts slot
#' @param count_type Required if `moo_counts` is a `multiOmicDataSet`: the type of counts to use -- must be a name in the counts slot (`moo@counts`).
#' @param sub_count_type Used if `moo_counts` is a `multiOmicDataSet` AND if `count_type` is a list, specify the sub count type within the list
#' @param ... arguments forwarded to method
#'
#' @returns ggplot object
#' @export
#'
#' @examples
#' # plot histogram for a counts slot in a multiOmicDataset Object
#' moo <- multiOmicDataSet(
#'   sample_metadata = nidap_sample_metadata,
#'   anno_dat = data.frame(),
#'   counts_lst = list("raw" = nidap_raw_counts)
#' )
#' p <- plot_histogram(moo, count_type = "raw")
#'
#' # customize the plot
#' plot_histogram(moo,
#'   count_type = "raw",
#'   group_colname = "GeneName", color_by_group = TRUE
#' )
#'
#' # plot histogram for a counts dataframe directly
#' counts_dat <- moo@counts$raw
#' plot_histogram(
#'   counts_dat,
#'   sample_metadata = nidap_sample_metadata,
#'   sample_id_colname = "Sample",
#'   feature_id_colname = "GeneName",
#'   label_colname = "Label"
#' )
#'
#' @details
#'
#' # Methods
#'
#' | link to docs  | class  |
#' |---|---|
#' | [plot_histogram_moo] | `multiOmicDataSet` |
#' | [plot_histogram_dat] | `data.frame`       |
#'
#' @name plot_histogram
#' @rdname plot_histogram
#' @family plotters
plot_histogram <- S7::new_generic("plot_histogram", dispatch_args = "moo_counts")

#' Plot histogram for multiOmicDataSet
#'
#' @param moo_counts counts dataframe or `multiOmicDataSet` containing `count_type` & `sub_count_type` in the counts slot
#' @param count_type Required if `moo_counts` is a `multiOmicDataSet`: the type of counts to use -- must be a name in the counts slot (`moo@counts`).
#' @param sub_count_type Used if `moo_counts` is a `multiOmicDataSet` AND if `count_type` is a list, specify the sub count type within the list
#' @param ... arguments forwarded to method: [plot_histogram_dat]
#' @examples
#' # plot histogram for a counts slot in a multiOmicDataset Object
#' moo <- multiOmicDataSet(
#'   sample_metadata = nidap_sample_metadata,
#'   anno_dat = data.frame(),
#'   counts_lst = list("raw" = nidap_raw_counts)
#' )
#' p <- plot_histogram(moo, count_type = "raw")
#'
#' # customize the plot
#' plot_histogram(moo,
#'   count_type = "raw",
#'   group_colname = "Gene", color_by_group = TRUE
#' )
#'
#' @name plot_histogram_moo
#' @method plot_histogram multiOmicDataSet
#' @seealso [plot_histogram] generic
#' @family plotters for multiOmicDataSets
S7::method(plot_histogram, multiOmicDataSet) <- function(moo_counts,
                                                         count_type,
                                                         sub_count_type = NULL,
                                                         ...) {
  counts_dat <- extract_counts(moo_counts, count_type, sub_count_type)
  plot_histogram(counts_dat, sample_metadata = moo_counts@sample_meta, ...)
}

#' Plot histogram for counts dataframe
#'
#' @param moo_counts counts dataframe (**required**)
#' @param sample_metadata sample metadata as a data frame or tibble (**required**)
#' @param sample_id_colname The column from the sample metadata containing the sample names. The names in this column must exactly match the names used as the sample column names of your input Counts Matrix. (Default: `NULL` - first column in the sample metadata will be used.)
#' @param feature_id_colname The column from the counts dataa containing the Feature IDs (Usually Gene or Protein ID). This is usually the first column of your input Counts Matrix. Only columns of Text type from your input Counts Matrix will be available to select for this parameter. (Default: `NULL` - first column in the counts matrix will be used.)
#' @param group_colname The column from the sample metadata containing the sample group information. This is usually a column showing to which experimental treatments each sample belongs (e.g. WildType, Knockout, Tumor, Normal, Before, After, etc.).
#' @param label_colname The column from the sample metadata containing the sample labels as you wish them to appear in the plots produced by this template. This can be the same Sample Names Column. However, you may desire different labels to display on your figure (e.g. shorter labels are sometimes preferred on plots). In that case, select the column with your preferred Labels here. The selected column should contain unique names for each sample. (Default: `NULL` -- `sample_id_colname` will be used.)
#' @param color_values vector of colors as hex values or names recognized by R
#' @param color_by_group Set to FALSE to label histogram by Sample Names, or set to TRUE to label histogram by the column you select in the "Group Column Used to Color Histogram" parameter (below). Default is FALSE.
#' @param set_min_max_for_x_axis whether to override the default for `ggplot2::xlim()` (default: `FALSE`)
#' @param minimum_for_x_axis value to override default `min` for `ggplot2::xlim()`
#' @param maximum_for_x_axis value to override default `max` for `ggplot2::xlim()`
#' @param x_axis_label text label for the x axis `ggplot2::xlab()`
#' @param y_axis_label text label for the y axis `ggplot2::ylab()`
#' @param legend_position passed to in `legend.position` `ggplot2::theme()`
#' @param legend_font_size passed to `ggplot2::element_text()` via `ggplot2::theme()`
#' @param number_of_legend_columns passed to `ncol` in `ggplot2::guide_legend()`
#' @param interactive_plots set to TRUE to make the plot interactive with `plotly`, allowing you to hover your mouse over a point or line to view sample information. The similarity heat map will not display if this toggle is set to TRUE. Default is FALSE.
#' @examples
#'
#' # plot histogram for a counts dataframe directly
#' counts_dat <- moo@counts$raw
#' plot_histogram(
#'   counts_dat,
#'   sample_metadata = nidap_sample_metadata,
#'   sample_id_colname = "Sample",
#'   feature_id_colname = "GeneName",
#'   label_colname = "Label"
#' )
#'
#' # customize the plot
#' plot_histogram(
#'   counts_dat,
#'   sample_metadata = nidap_sample_metadata,
#'   sample_id_colname = "Sample",
#'   feature_id_colname = "GeneName",
#'   group_colname = "Group",
#'   color_by_group = TRUE
#' )
#'
#' @seealso [plot_histogram] generic
#' @name plot_histogram_dat
#' @method plot_histogram data.frame
#' @family plotters for counts dataframes
S7::method(plot_histogram, S7::class_data.frame) <- function(moo_counts,
                                                             sample_metadata,
                                                             sample_id_colname = NULL,
                                                             feature_id_colname = NULL,
                                                             group_colname = "Group",
                                                             label_colname = "Label",
                                                             color_values = c(
                                                               "#5954d6",
                                                               "#e1562c",
                                                               "#b80058",
                                                               "#00c6f8",
                                                               "#d163e6",
                                                               "#00a76c",
                                                               "#ff9287",
                                                               "#008cf9",
                                                               "#006e00",
                                                               "#796880",
                                                               "#FFA500",
                                                               "#878500"
                                                             ),
                                                             color_by_group = FALSE,
                                                             set_min_max_for_x_axis = FALSE,
                                                             minimum_for_x_axis = -1,
                                                             maximum_for_x_axis = 1,
                                                             x_axis_label = "Counts",
                                                             y_axis_label = "Density",
                                                             legend_position = "top",
                                                             legend_font_size = 10,
                                                             number_of_legend_columns = 6,
                                                             interactive_plots = FALSE) {
  count <- NULL
  counts_dat <- moo_counts
  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(counts_dat)[1]
  }

  df_long <- counts_dat %>%
    tidyr::pivot_longer(
      -tidyselect::all_of(feature_id_colname),
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
      ggplot2::ggplot(ggplot2::aes(
        x = count,
        group = !!rlang::sym(sample_id_colname)
      )) +
      ggplot2::geom_density(ggplot2::aes(colour = !!rlang::sym(group_colname)), linewidth = 1)
  } else {
    n <- df_long %>%
      dplyr::pull(sample_id_colname) %>%
      unique() %>%
      length()

    hist_plot <- df_long %>%
      ggplot2::ggplot(ggplot2::aes(
        x = count,
        group = !!rlang::sym(sample_id_colname)
      )) +
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
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0
      ),
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
