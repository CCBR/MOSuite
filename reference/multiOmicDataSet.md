# multiOmicDataSet class

multiOmicDataSet class

## Usage

``` r
multiOmicDataSet(sample_metadata, anno_dat, counts_lst, analyses_lst = list())
```

## Arguments

- sample_metadata:

  sample metadata as a data frame or tibble. The first column is assumed
  to contain the sample IDs which must correspond to column names in the
  raw counts.

- anno_dat:

  data frame of feature annotations, such as gene symbols or any other
  information about the features in `counts_lst`.

- counts_lst:

  named list of data frames containing counts, e.g. expected feature
  counts from RSEM. Each data frame is expected to contain a
  `feature_id` column as the first column, and all remaining columns are
  sample IDs in the `sample_meta`.

- analyses_lst:

  named list of analysis results, e.g. DESeq results object

## See also

Other moo constructors:
[`create_multiOmicDataSet_from_dataframes()`](https://ccbr.github.io/MOSuite/reference/create_multiOmicDataSet_from_dataframes.md),
[`create_multiOmicDataSet_from_files()`](https://ccbr.github.io/MOSuite/reference/create_multiOmicDataSet_from_files.md)
