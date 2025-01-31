% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/histogram.R
\name{plot_histogram}
\alias{plot_histogram}
\title{Plot histogram}
\usage{
plot_histogram(
  log_counts,
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
  make_plots_interactive = FALSE
)
}
\arguments{
\item{log_counts}{log-transformed filtered counts}

\item{sample_metadata}{sample metadata as a data frame or tibble.}

\item{feature_id_colname}{The column from the counts dataa containing the Feature IDs (Usually Gene or Protein ID). This is usually the first column of your input Counts Matrix. Only columns of Text type from your input Counts Matrix will be available to select for this parameter. (Default: \code{NULL} - first column in the counts matrix will be used.)}

\item{group_colname}{The column from the sample metadata containing the sample group information. This is usually a column showing to which experimental treatments each sample belongs (e.g. WildType, Knockout, Tumor, Normal, Before, After, etc.).}

\item{label_colname}{The column from the sample metadata containing the sample labels as you wish them to appear in the plots produced by this template. This can be the same Sample Names Column. However, you may desire different labels to display on your figure (e.g. shorter labels are sometimes preferred on plots). In that case, select the column with your preferred Labels here. The selected column should contain unique names for each sample. (Default: \code{NULL} -- \code{sample_id_colname} will be used.)}

\item{color_values}{vector of colors as hex values or names recognized by R}

\item{color_histogram_by_group}{Set to FALSE to label histogram by Sample Names, or set to TRUE to label histogram by the column you select in the "Group Column Used to Color Histogram" parameter (below). Default is FALSE.}

\item{set_min_max_for_x_axis_for_histogram}{whether to set min/max value for histogram x-axis}

\item{minimum_for_x_axis_for_histogram}{x-axis minimum for histogram plot}

\item{maximum_for_x_axis_for_histogram}{x-axis maximum for histogram plot}

\item{legend_position_for_histogram}{legend position for the histogram plot. consider setting to 'none' for a large number of samples.}

\item{legend_font_size_for_histogram}{legend font size for the histogram plot}

\item{number_of_histogram_legend_columns}{number of columns for the histogram legend}

\item{x_axis_label}{text label for the x axis}

\item{y_axis_label}{text label for the y axis}

\item{make_plots_interactive}{set to TRUE to make PCA and Histogram plots interactive with \code{plotly}, allowing you to hover your mouse over a point or line to view sample information. The similarity heat map will not display if this toggle is set to TRUE. Default is FALSE.}
}
\value{
ggplot object
}
\description{
Plot histogram
}
