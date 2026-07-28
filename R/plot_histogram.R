#' Plot histogram
#'
#' @param moo_counts counts dataframe or `multiOmicDataSet` containing `count_type` & `sub_count_type` in the counts
#'   slot
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
#'   group_colname = "Group", color_by_group = TRUE
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
#' @seealso
#' - [plot_histogram.multiOmicDataSet()]
#' - [plot_histogram.data.frame()]
#'
#' @family plotters
#' @keywords plotters
#' @family moo methods
plot_histogram <- S7::new_generic(
  "plot_histogram",
  dispatch_args = "moo_counts"
)

#' Plot histogram for multiOmicDataSet
#'
#' @rdname plot_histogram.multiOmicDataSet
#' @aliases plot_histogram.multiOmicDataSet
#' @usage NULL
#'
#' @param count_type Required if `moo_counts` is a `multiOmicDataSet`: the type of counts to use -- must be a name in
#'   the counts slot (`moo@counts`).
#' @param sub_count_type Used if `moo_counts` is a `multiOmicDataSet` AND if `count_type` is a list, specify the sub
#'   count type within the list
#' @param group_colname The column from the sample metadata containing the sample group information. This is usually a
#'   column showing to which experimental treatments each sample belongs (e.g. WildType, Knockout, Tumor, Normal,
#'   Before, After, etc.).
#' @param color_values vector of colors as hex values or names recognized by R. Unnamed colors are assigned by factor
#'   level order when the grouping column is a factor; otherwise, they follow the order in which groups first appear in
#'   the metadata column. Defaults to `NULL`; when `NULL`, `mosuite_palette` is used for `data.frame` dispatch and
#'   stored colors are used for `multiOmicDataSet` dispatch.
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
#'   group_colname = "Group", color_by_group = TRUE
#' )
#'
#' @seealso [plot_histogram()] generic
#' @family plotters for multiOmicDataSets
S7::method(plot_histogram, multiOmicDataSet) <- function(
  moo_counts,
  count_type,
  sub_count_type = NULL,
  group_colname = "Group",
  color_values = NULL,
  ...
) {
  counts_dat <- extract_counts(moo_counts, count_type, sub_count_type)
  color_values <- color_values %||% moo_counts@analyses$colors[[group_colname]]
  return(plot_histogram(
    counts_dat,
    sample_metadata = moo_counts@sample_meta,
    group_colname = group_colname,
    color_values = color_values,
    ...
  ))
}

build_histogram_hover_text <- function(
  histogram_data,
  sample_id_colname,
  group_colname = NULL
) {
  return(format_hover_text(
    histogram_data,
    primary_colname = sample_id_colname,
    secondary_colname = group_colname,
    missing_col_context = "histogram",
    require_secondary = FALSE
  ))
}

#' Plot histogram for counts dataframe
#'
#' @rdname plot_histogram.data.frame
#' @aliases plot_histogram.data.frame
#' @usage NULL
#'
#' @param sample_metadata sample metadata as a data frame or tibble (**required**)
#' @param sample_id_colname The column from the sample metadata containing the sample names. The names in this column
#'   must exactly match the names used as the sample column names of your input Counts Matrix. (Default: `NULL` - first
#'   column in the sample metadata will be used.)
#' @param feature_id_colname The column from the counts dataa containing the Feature IDs (Usually Gene or Protein ID).
#'   This is usually the first column of your input Counts Matrix. Only columns of Text type from your input Counts
#'   Matrix will be available to select for this parameter. (Default: `NULL` - first column in the counts matrix will be
#'   used.)
#' @param group_colname The column from the sample metadata containing the sample group information. This is usually a
#'   column showing to which experimental treatments each sample belongs (e.g. WildType, Knockout, Tumor, Normal,
#'   Before, After, etc.).
#' @param label_colname The column from the sample metadata containing the sample labels as you wish them to appear in
#'   the plots produced by this template. This can be the same Sample Names Column. However, you may desire different
#'   labels to display on your figure (e.g. shorter labels are sometimes preferred on plots). In that case, select the
#'   column with your preferred Labels here. The selected column should contain unique names for each sample. (Default:
#'   `NULL` -- `sample_id_colname` will be used.)
#' @param color_values vector of colors as hex values or names recognized by R. Unnamed colors are assigned by factor
#'   level order when the grouping column is a factor; otherwise, they follow the order in which groups first appear in
#'   the metadata column. Defaults to `NULL`; when `NULL`, `mosuite_palette` is used.
#' @param color_by_group Set to FALSE to label histogram by Sample Names, or set to TRUE to label histogram by the
#'   column you select in the "Group Column Used to Color Histogram" parameter (below). Default is FALSE.
#' @param set_min_max_for_x_axis whether to override the default for `ggplot2::xlim()` (default: `FALSE`)
#' @param minimum_for_x_axis value to override default `min` for `ggplot2::xlim()`
#' @param maximum_for_x_axis value to override default `max` for `ggplot2::xlim()`
#' @param x_axis_label text label for the x axis `ggplot2::xlab()`
#' @param y_axis_label text label for the y axis `ggplot2::ylab()`
#' @param legend_position passed to in `legend.position` `ggplot2::theme()`
#' @param legend_font_size passed to `ggplot2::element_text()` via `ggplot2::theme()`. If `NULL`, the size is scaled
#'   automatically based on the number and length of legend labels.
#' @param number_of_legend_columns passed to `ncol` in `ggplot2::guide_legend()`
#' @param interactive_plots set to TRUE to make the plot interactive with `plotly`, allowing you to hover your mouse
#'   over a point or line to view sample information. The similarity heat map will not display if this toggle is set to
#'   TRUE. Default is FALSE.
#' @param return_ggplot If `TRUE`, return the ggplot object prepared for interactive hover text before converting it to
#'   plotly. Used when callers need to add more ggplot layers first. Default is `FALSE`.
#' @param use_log2_x_axis If `TRUE`, add a display-only pseudocount to plotted values and use a log2 x-axis. Default is
#'   `FALSE`.
#' @param ... additional arguments (ignored; accepted for compatibility with the moo dispatch)
#' @examples
#'
#' # plot histogram for a counts dataframe directly
#' plot_histogram(
#'   nidap_clean_raw_counts,
#'   sample_metadata = nidap_sample_metadata,
#'   sample_id_colname = "Sample",
#'   feature_id_colname = "Gene",
#'   label_colname = "Label"
#' )
#'
#' # customize the plot
#' plot_histogram(
#'   nidap_clean_raw_counts,
#'   sample_metadata = nidap_sample_metadata,
#'   sample_id_colname = "Sample",
#'   feature_id_colname = "Gene",
#'   group_colname = "Group",
#'   color_by_group = TRUE
#' )
#'
#' @seealso [plot_histogram()] generic
#'
#' @family plotters for counts dataframes
S7::method(plot_histogram, S7::class_data.frame) <- function(
  moo_counts,
  sample_metadata,
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  group_colname = "Group",
  label_colname = "Label",
  color_values = NULL,
  color_by_group = FALSE,
  set_min_max_for_x_axis = FALSE,
  minimum_for_x_axis = -1,
  maximum_for_x_axis = 1,
  x_axis_label = "Counts",
  y_axis_label = "Density",
  legend_position = "top",
  legend_font_size = NULL,
  number_of_legend_columns = 6,
  interactive_plots = FALSE,
  return_ggplot = FALSE,
  use_log2_x_axis = FALSE,
  ...
) {
  count <- NULL
  log2_axis_pseudocount <- 0.5
  color_values <- color_values %||% mosuite_palette
  counts_dat <- moo_counts
  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(counts_dat)[1]
  }

  df_long <- counts_dat |>
    tidyr::pivot_longer(
      -tidyselect::all_of(feature_id_colname),
      names_to = sample_id_colname,
      values_to = "count"
    ) |>
    dplyr::left_join(sample_metadata, by = sample_id_colname)

  # For log2 histogram axes, add a display-only pseudocount after reshaping.
  # The original count table is not modified; this only keeps zeros finite for ggplot's log transform.
  if (isTRUE(use_log2_x_axis)) {
    df_long <- df_long |>
      dplyr::mutate(count = count + log2_axis_pseudocount)
  }

  # Match user-supplied x-axis limits to the plotted scale. On the log2 path,
  # limits need the same pseudocount offset as the displayed data.
  if (set_min_max_for_x_axis == TRUE) {
    if (isTRUE(use_log2_x_axis)) {
      xmin <- minimum_for_x_axis + log2_axis_pseudocount
      xmax <- maximum_for_x_axis + log2_axis_pseudocount
    } else {
      xmin <- minimum_for_x_axis
      xmax <- maximum_for_x_axis
    }
  } else {
    xmin <- min(df_long |> dplyr::pull(count))
    xmax <- max(df_long |> dplyr::pull(count))
  }
  # Guard the log2 axis against non-finite or nonpositive limits after the offset.
  if (isTRUE(use_log2_x_axis)) {
    if (set_min_max_for_x_axis == TRUE) {
      xmin <- max(xmin, log2_axis_pseudocount)
    } else {
      # Automatic log2 histograms should display a raw minimum of 0.
      xmin <- log2_axis_pseudocount
    }
    if (!is.finite(xmax) || xmax <= xmin) {
      xmax <- max(xmin, log2_axis_pseudocount)
    }
  }

  if (color_by_group == TRUE) {
    df_long <- df_long |>
      dplyr::filter(!is.na(!!rlang::sym(group_colname)))
    color_values <- resolve_plot_colors(df_long, group_colname, color_values)
    df_long <- df_long |>
      dplyr::mutate(
        !!rlang::sym(group_colname) := as.character(!!rlang::sym(group_colname))
      )
    if (isTRUE(interactive_plots)) {
      df_long$histogram_hover_text <- build_histogram_hover_text(
        df_long,
        sample_id_colname,
        group_colname
      )
    }
    histogram_mapping <- ggplot2::aes(
      x = count,
      group = !!rlang::sym(sample_id_colname)
    )
    if (isTRUE(interactive_plots)) {
      histogram_mapping <- ggplot2::aes(
        x = count,
        group = !!rlang::sym(sample_id_colname),
        text = histogram_hover_text
      )
    }
    # The problem here is that static histograms should keep density curves in
    # the plot while showing line-style legend keys instead of box-like keys.
    # We build geom_density() args as a list so static output can add
    # key_glyph, but the interactive ggplotly() path can skip that tweak.
    # Passing the legend-key change into ggplotly() made interactive density
    # traces misbehave, so we apply it only for non-interactive plots.
    density_layer_args <- list(
      mapping = ggplot2::aes(colour = !!rlang::sym(group_colname)),
      linewidth = 1
    )
    if (!isTRUE(interactive_plots)) {
      density_layer_args$key_glyph <- ggplot2::draw_key_path
    }

    # plot Density
    hist_plot <- df_long |>
      ggplot2::ggplot(histogram_mapping) +
      do.call(ggplot2::geom_density, density_layer_args)
  } else {
    color_values <- resolve_plot_colors(
      df_long,
      sample_id_colname,
      color_values
    )
    df_long <- df_long |>
      dplyr::mutate(
        !!rlang::sym(sample_id_colname) := as.character(
          !!rlang::sym(sample_id_colname)
        )
      )
    if (isTRUE(interactive_plots)) {
      df_long$histogram_hover_text <- build_histogram_hover_text(
        df_long,
        sample_id_colname,
        group_colname
      )
    }
    histogram_mapping <- ggplot2::aes(
      x = count,
      group = !!rlang::sym(sample_id_colname)
    )
    if (isTRUE(interactive_plots)) {
      histogram_mapping <- ggplot2::aes(
        x = count,
        group = !!rlang::sym(sample_id_colname),
        text = histogram_hover_text
      )
    }
    # Use the same strategy for sample-colored histograms: solve the static
    # legend-key problem without changing the interactive density conversion.
    # Static output gets line-style legend keys, while interactive output
    # avoids the key_glyph change that breaks ggplotly().
    density_layer_args <- list(
      mapping = ggplot2::aes(colour = !!rlang::sym(sample_id_colname)),
      linewidth = 1
    )
    if (!isTRUE(interactive_plots)) {
      density_layer_args$key_glyph <- ggplot2::draw_key_path
    }

    hist_plot <- df_long |>
      ggplot2::ggplot(histogram_mapping) +
      do.call(ggplot2::geom_density, density_layer_args)
  }

  legend_font_size <- get_legend_text_size(
    names(color_values),
    legend_font_size
  )

  # Keep the plain histogram default unchanged. When requested, use ggplot's log2
  # scale and label ticks back on the original count scale by subtracting the offset.
  x_axis_scale <- if (isTRUE(use_log2_x_axis)) {
    ggplot2::scale_x_continuous(
      transform = "log2",
      limits = c(xmin, xmax),
      breaks = function(limits) {
        # Breaks are chosen on the original count scale, then shifted by the
        # display-only pseudocount so they align with the plotted/log2 values.
        raw_upper_limit <- max(limits - log2_axis_pseudocount, na.rm = TRUE)
        if (!is.finite(raw_upper_limit) || raw_upper_limit <= 0) {
          return(log2_axis_pseudocount)
        }
        max_power <- floor(log2(raw_upper_limit))
        if (max_power < 0) {
          raw_breaks <- c(0, raw_upper_limit)
        } else {
          # Use powers of two, thinning to about eight labels for wide ranges,
          # and always include the highest power so the high end is labeled.
          step <- max(1, ceiling((max_power + 1) / 8))
          exponents <- unique(c(seq(0, max_power, by = step), max_power))
          raw_breaks <- c(0, 2^exponents)
        }
        raw_breaks <- raw_breaks[
          raw_breaks >= 0 & raw_breaks <= raw_upper_limit
        ]
        return(unique(raw_breaks + log2_axis_pseudocount))
      },
      labels = function(x) {
        return(scales::label_number(big.mark = "")(x - log2_axis_pseudocount))
      },
      name = x_axis_label
    )
  } else {
    ggplot2::xlim(xmin, xmax)
  }

  hist_plot <- hist_plot +
    ggplot2::xlab(x_axis_label) +
    ggplot2::ylab(y_axis_label) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = legend_position,
      legend.key = ggplot2::element_blank(),
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
    x_axis_scale +
    # scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),n)) +
    ggplot2::scale_colour_manual(values = color_values)

  if (isTRUE(interactive_plots)) {
    if (isTRUE(return_ggplot)) {
      return(hist_plot)
    }
    hist_plot <- hist_plot |>
      plotly::ggplotly(tooltip = "text")
  } else {
    hist_plot <- add_colour_legend_layout(
      hist_plot,
      labels = names(color_values),
      legend_position = legend_position,
      ncol = number_of_legend_columns,
      legend_text_size = legend_font_size,
      guide_override_aes = list(
        linetype = 1,
        linewidth = 2,
        shape = NA,
        fill = NA
      )
    )
  }
  return(hist_plot)
}
