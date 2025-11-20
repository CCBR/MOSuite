# Calculate CPM on a data frame

Calculate CPM on a data frame

## Usage

``` r
calc_cpm_df(dat, feature_id_colname = "gene_id", ...)
```

## Arguments

- dat:

  data frame of counts with a gene column

- feature_id_colname:

  name of the column in `counts_dat` that contains feature/gene IDs.
  (Default: `NULL` - first column in the count data will be used.)

- ...:

  additional arguments to pass to edger::cpm()

## Value

cpm-transformed counts as a data frame
