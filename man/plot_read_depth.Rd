% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_read_depth.R
\name{plot_read_depth}
\alias{plot_read_depth}
\title{Plot read depth as a bar plot}
\usage{
plot_read_depth(moo_counts, ...)
}
\arguments{
\item{moo_counts}{counts dataframe or \code{multiOmicDataSet} containing \code{count_type} & \code{sub_count_type} in the counts slot}

\item{...}{arguments forwarded to method}
}
\value{
ggplot barplot
}
\description{
The first argument can be a \code{multiOmicDataset} object (\code{moo}) or a \code{data.frame} containing counts.
For a \code{moo}, choose which counts slot to use with \code{count_type} & (optionally) \code{sub_count_type}.
}
\section{Methods}{
\tabular{ll}{
   link to docs \tab class \cr
   \link{plot_read_depth_moo} \tab \code{multiOmicDataSet} \cr
   \link{plot_read_depth_dat} \tab \code{data.frame} \cr
}
}

\examples{
# multiOmicDataSet
moo <- multiOmicDataSet(
  sample_metadata = nidap_sample_metadata,
  anno_dat = data.frame(),
  counts_lst = list(
    "raw" = nidap_raw_counts,
    "clean" = nidap_clean_raw_counts
  )
)

plot_read_depth(moo, count_type = "clean")

# dataframe
plot_read_depth(nidap_clean_raw_counts)

}
\seealso{
Other plotters: 
\code{\link{plot_corr_heatmap}()},
\code{\link{plot_expr_heatmap}()},
\code{\link{plot_histogram}()},
\code{\link{plot_pca}()},
\code{\link{print_or_save_plot}()}

Other moo methods: 
\code{\link{batch_correct_counts}()},
\code{\link{clean_raw_counts}()},
\code{\link{diff_counts}()},
\code{\link{filter_counts}()},
\code{\link{filter_diff}()},
\code{\link{normalize_counts}()},
\code{\link{plot_corr_heatmap}()},
\code{\link{plot_expr_heatmap}()},
\code{\link{plot_histogram}()},
\code{\link{plot_pca}()},
\code{\link{run_deseq2}()},
\code{\link{set_color_pal}()}
}
\concept{moo methods}
\concept{plotters}
\keyword{plotters}
