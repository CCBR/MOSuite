% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{print_or_save_plot}
\alias{print_or_save_plot}
\title{Print and/or save a ggplot}
\usage{
print_or_save_plot(
  plot_obj,
  filename,
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_dir = options::opt("plots_dir"),
  graphics_device = grDevices::png,
  ...
)
}
\arguments{
\item{plot_obj}{plot object (e.g. ggplot, ComplexHeatmap...)}

\item{filename}{name of the output file. will be joined with the \code{plots_dir} option.}

\item{print_plots}{Whether to print plots during analysis (Defaults to \code{FALSE}, overwritable using option 'moo_print_plots' or environment variable 'MOO_PRINT_PLOTS')}

\item{save_plots}{Whether to save plots to files during analysis (Defaults to \code{FALSE}, overwritable using option 'moo_save_plots' or environment variable 'MOO_SAVE_PLOTS')}

\item{plots_dir}{Path where plots are saved when \code{moo_save_plots} is \code{TRUE} (Defaults to \code{"figures/"}, overwritable using option 'moo_plots_dir' or environment variable 'MOO_PLOTS_DIR')}

\item{graphics_device}{Default: \code{grDevice::png()}. Only used if the plot is not a ggplot.}

\item{...}{arguments forwarded to \code{ggplot2::ggsave()}}
}
\value{
invisibly returns the path where the plot image was saved to the disk
}
\description{
If \code{save_plots} is \code{TRUE}, the plot will be saved as an image to the path at
\code{file.path(plots_dir, filename)}.
If \code{plot_obj} is a ggplot, \code{ggplot2::ggsave()} is used to save the image.
Otherwise, \code{graphics_device} is used (\code{grDevice::png()} by default).
}
\seealso{
Other plotters: 
\code{\link{plot_corr_heatmap}()},
\code{\link{plot_expr_heatmap}()},
\code{\link{plot_histogram}()},
\code{\link{plot_pca}()},
\code{\link{plot_read_depth}()}
}
\concept{plotters}
\keyword{plotters}
