% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean.R
\name{check_sample_names}
\alias{check_sample_names}
\title{check sample names in metadata & counts data}
\usage{
check_sample_names(counts, metadata, sample_names_column)
}
\arguments{
\item{counts}{dataframe containing counts data}

\item{metadata}{dataframe containing sample metadata}

\item{sample_names_column}{The column from your input Sample Metadata table containing the sample names. The names in this column must exactly match the names used as the sample column names of your input Counts Matrix. Only columns of Text type from your input Sample Metadata table will be available to select for this parameter.}
}
\description{
check sample names in metadata & counts data
}
