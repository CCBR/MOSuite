% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0_mo-class.R
\name{create_multiOmicDataSet_from_files}
\alias{create_multiOmicDataSet_from_files}
\title{Construct a multiOmicDataSet object from tsv files.}
\usage{
create_multiOmicDataSet_from_files(
  sample_meta_filepath,
  feature_counts_filepath,
  count_type = "raw",
  sample_id_colname = NULL,
  feature_id_colname = NULL
)
}
\arguments{
\item{sample_meta_filepath}{path to tsv file with sample IDs and metadata for differential analysis.}

\item{feature_counts_filepath}{path to tsv file of expected feature counts (e.g. gene counts from RSEM).}

\item{count_type}{type to assign the values of \code{counts_dat} to in the \code{counts} slot}

\item{sample_id_colname}{name of the column in \code{sample_metadata} that contains the sample IDs. (Default: \code{NULL} - first column in the sample metadata will be used.)}

\item{feature_id_colname}{name of the column in \code{counts_dat} that contains feature/gene IDs. (Default: \code{NULL} - first column in the count data will be used.)}
}
\value{
\link{multiOmicDataSet} object
}
\description{
Construct a multiOmicDataSet object from tsv files.
}
\examples{
moo <- create_multiOmicDataSet_from_files(
  sample_meta_filepath = system.file("extdata",
    "sample_metadata.tsv.gz",
    package = "MOSuite"
  ),
  feature_counts_filepath = system.file("extdata",
    "RSEM.genes.expected_count.all_samples.txt.gz",
    package = "MOSuite"
  )
)
moo@counts$raw \%>\% head()
moo@sample_meta

}
\seealso{
Other moo constructors: 
\code{\link{create_multiOmicDataSet_from_dataframes}()},
\code{\link{multiOmicDataSet}()}
}
\concept{moo constructors}
