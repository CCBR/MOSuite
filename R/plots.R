#' Print and/or save a ggplot
#'
#' If `save_plots` is `TRUE`, the plot will be saved as an image to the path at
#' `file.path(plots_dir, filename)`.
#' If `plot_obj` is a ggplot, `ggplot2::ggsave()` is used to save the image.
#' Otherwise, `graphics_device` is used (`grDevices::png()` by default).
#'
#' @inheritParams option_params
#' @param plot_obj plot object (e.g. ggplot, ComplexHeatmap...)
#' @param filename name of the output file. will be joined with the `plots_dir` option.
#' @param graphics_device Default: `grDevices::png()`. Only used if the plot is not a ggplot.
#' @param caption optional caption text to add to the plot. For ggplot objects, this is
#'   added via `ggplot2::labs(caption = caption)`. For `ComplexHeatmap` objects, the
#'   caption is rendered at the bottom of the graphics device using `grid::grid.text()`.
#' @param ... arguments forwarded to `ggplot2::ggsave()`
#'
#' @return invisibly returns the path where the plot image was saved to the disk
#' @export
#' @family plotters
#' @keywords plotters
print_or_save_plot <- function(
  plot_obj,
  filename,
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_dir = options::opt("plots_dir"),
  graphics_device = grDevices::png,
  caption = NULL,
  ...
) {
  draw_heatmap_with_caption <- function(hm) {
    ComplexHeatmap::draw(hm)
    if (!is.null(caption)) {
      grid::grid.text(
        caption,
        x = grid::unit(0.5, "npc"),
        y = grid::unit(2, "mm"),
        just = "bottom",
        gp = grid::gpar(fontsize = 9, col = "grey40")
      )
    }
  }
  if (!is.null(caption) && inherits(plot_obj, "ggplot")) {
    plot_obj <- plot_obj + ggplot2::labs(caption = caption)
  }
  if (isTRUE(print_plots)) {
    if (inherits(plot_obj, c("Heatmap", "HeatmapList"))) {
      draw_heatmap_with_caption(plot_obj)
    } else {
      print(plot_obj)
    }
  }
  if (isTRUE(save_plots)) {
    # create output directory if it doesn't exist
    if (!is.null(plots_dir) && nchar(plots_dir) > 0) {
      filename <- file.path(plots_dir, filename)
    }
    outdir <- dirname(filename)
    if (!dir.exists(outdir)) {
      dir.create(outdir, recursive = TRUE)
    }

    # select saving methods depending on plot object class
    if (inherits(plot_obj, "ggplot")) {
      ggplot2::ggsave(filename = filename, plot = plot_obj, ...)
    } else if (inherits(plot_obj, "htmlwidget")) {
      htmlwidgets::saveWidget(plot_obj, filename, ...)
    } else if (inherits(plot_obj, c("Heatmap", "HeatmapList"))) {
      graphics_device(file = filename)
      on.exit(grDevices::dev.off(), add = TRUE)
      draw_heatmap_with_caption(plot_obj)
    } else {
      graphics_device(file = filename)
      on.exit(grDevices::dev.off(), add = TRUE)
      plot(plot_obj)
    }
  }
  return(invisible(filename))
}

#' Compute a wrapped colour legend column count
#'
#' Computes a conservative number of legend columns for horizontal ggplot colour
#' legends. Top and bottom legends are wrapped based on the number of labels and
#' the longest label length. Other legend positions return `NULL` so their
#' existing ggplot layout is preserved.
#'
#' @param labels Character vector of legend labels.
#' @param legend_position Legend position passed to `ggplot2::theme()`.
#' @param ncol Optional maximum number of legend columns.
#' @param legend_text_size Legend text size used to scale the horizontal space
#'   estimate. Larger legend text uses fewer columns.
#' @param max_label_characters_per_row Approximate total label characters to fit
#'   on one horizontal legend row.
#'
#' @return Integer column count for top/bottom legends, or `NULL` when no
#'   wrapping should be applied.
#' @keywords internal
#' @noRd
get_legend_column_count <- function(
  labels,
  legend_position = "top",
  ncol = NULL,
  legend_text_size = 10,
  max_label_characters_per_row = 45
) {
  if (!legend_position %in% c("top", "bottom")) {
    return(NULL)
  }

  labels <- stats::na.omit(as.character(labels))
  if (length(labels) == 0) {
    return(NULL)
  }

  max_label_length <- max(nchar(labels), 1)
  text_size_multiplier <- legend_text_size / 10
  columns_by_label_length <- max(
    1,
    floor(
      max_label_characters_per_row / (max_label_length * text_size_multiplier)
    )
  )
  legend_columns <- min(length(labels), columns_by_label_length)
  if (!is.null(ncol)) {
    legend_columns <- min(ncol, legend_columns)
  }
  legend_columns
}

#' Compute colour legend text size
#'
#' Computes legend text size from legend labels. Short legends keep the larger
#' default text used by simple group legends, while longer or denser legends are
#' scaled down.
#'
#' @param labels Character vector of legend labels.
#' @param legend_text_size Optional explicit legend text size. When supplied,
#'   this value is returned unchanged.
#' @param min_legend_text_size Smallest automatically selected legend text size.
#' @param max_legend_text_size Largest automatically selected legend text size.
#'
#' @return Numeric legend text size.
#' @keywords internal
#' @noRd
get_legend_text_size <- function(
  labels,
  legend_text_size = NULL,
  min_legend_text_size = 8,
  max_legend_text_size = 18
) {
  if (!is.null(legend_text_size)) {
    return(legend_text_size)
  }

  labels <- stats::na.omit(as.character(labels))
  if (length(labels) == 0) {
    return(max_legend_text_size)
  }

  label_pressure <- max(
    max(nchar(labels), 1) / 5,
    length(labels) / 4,
    1
  )

  max(min_legend_text_size, floor(max_legend_text_size / sqrt(label_pressure)))
}

#' Add wrapped colour legend layout to a ggplot
#'
#' Applies a colour guide with a wrapped column count for top and bottom legends.
#' Left, right, and hidden legends are returned unchanged.
#'
#' @param plot A `ggplot2` plot object.
#' @inheritParams get_legend_column_count
#'
#' @return A `ggplot2` plot object with colour legend layout applied when needed.
#' @keywords internal
#' @noRd
add_colour_legend_layout <- function(
  plot,
  labels,
  legend_position = "top",
  ncol = NULL,
  legend_text_size = 10,
  max_label_characters_per_row = 45
) {
  legend_columns <- get_legend_column_count(
    labels = labels,
    legend_position = legend_position,
    ncol = ncol,
    legend_text_size = legend_text_size,
    max_label_characters_per_row = max_label_characters_per_row
  )

  if (is.null(legend_columns)) {
    return(plot)
  }

  plot +
    ggplot2::guides(
      colour = ggplot2::guide_legend(ncol = legend_columns, byrow = TRUE)
    ) +
    ggplot2::theme(
      legend.box = "vertical"
    )
}
