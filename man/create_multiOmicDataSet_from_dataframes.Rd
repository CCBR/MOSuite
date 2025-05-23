% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0_mo-class.R
\name{create_multiOmicDataSet_from_dataframes}
\alias{create_multiOmicDataSet_from_dataframes}
\title{Construct a multiOmicDataSet object from data frames}
\usage{
create_multiOmicDataSet_from_dataframes(
  sample_metadata,
  counts_dat,
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  count_type = "raw"
)
}
\arguments{
\item{sample_metadata}{sample metadata as a data frame or tibble.
The first column is assumed to contain the sample IDs which must correspond to column names in the raw counts.}

\item{counts_dat}{data frame of feature counts (e.g. expected feature counts from RSEM).}

\item{sample_id_colname}{name of the column in \code{sample_metadata} that contains the sample IDs. (Default: \code{NULL} - first column in the sample metadata will be used.)}

\item{feature_id_colname}{name of the column in \code{counts_dat} that contains feature/gene IDs. (Default: \code{NULL} - first column in the count data will be used.)}

\item{count_type}{type to assign the values of \code{counts_dat} to in the \code{counts} slot}
}
\value{
\link{multiOmicDataSet} object
}
\description{
Construct a multiOmicDataSet object from data frames
}
\examples{
sample_meta <- data.frame(
  sample_id = c("KO_S3", "KO_S4", "WT_S1", "WT_S2"),
  condition = factor(
    c("knockout", "knockout", "wildtype", "wildtype"),
    levels = c("wildtype", "knockout")
  )
)
moo <- create_multiOmicDataSet_from_dataframes(sample_meta, gene_counts)
head(moo@sample_meta)
head(moo@counts$raw)
head(moo@annotation)

}
\seealso{
Other moo constructors: 
\code{\link{create_multiOmicDataSet_from_files}()},
\code{\link{multiOmicDataSet}()}
}
\concept{moo constructors}
