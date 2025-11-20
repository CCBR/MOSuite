# Convert a data frame of gene counts to a matrix

Convert a data frame of gene counts to a matrix

## Usage

``` r
counts_dat_to_matrix(counts_tbl, feature_id_colname = NULL)
```

## Arguments

- counts_tbl:

  expected feature counts as a dataframe or tibble, with all columns
  except `feature_id_colname`

- feature_id_colname:

  name of the column in `counts_dat` that contains feature/gene IDs.
  (Default: `NULL` - first column in the count data will be used.)

## Value

matrix of gene counts with rows as gene IDs

## Examples

``` r
if (FALSE) { # \dontrun{
counts_dat_to_matrix(head(gene_counts))
} # }
```
