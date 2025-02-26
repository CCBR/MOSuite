#' Print and/or save a ggplot
#'
#' If `save_plots` is `TRUE`, the plot will be saved as an image to the path at
#' `file.path(plots_dir, filename)`.
#' If `plot_obj` is a ggplot, `ggplot2::ggsave()` is used to save the image.
#' Otherwise, `graphics_device` is used (`grDevice::png()` by default).
#'
#' @inheritParams option_params
#' @param plot_obj plot object (e.g. ggplot, ComplexHeatmap...)
#' @param filename name of the output file. will be joined with the `plots_dir` option.
#' @param graphics_device Default: `grDevice::png()`. Only used if the plot is not a ggplot.
#' @param ... arguments forwarded to `ggplot2::ggsave()`
#'
#' @return invisibly returns the path where the plot image was saved to the disk
#' @export
#' @family plotters
print_or_save_plot <- function(plot_obj,
                               filename,
                               print_plots = options::opt("print_plots"),
                               save_plots = options::opt("save_plots"),
                               plots_dir = options::opt("plots_dir"),
                               graphics_device = grDevices::png,
                               ...) {
  if (isTRUE(print_plots)) {
    print(plot_obj)
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
    } else {
      graphics_device(file = filename)
      plot(plot_obj)
      grDevices::dev.off()
    }
  }
  return(invisible(filename))
}
