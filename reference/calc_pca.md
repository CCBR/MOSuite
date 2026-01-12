# Perform principal components analysis

Perform principal components analysis

## Usage

``` r
calc_pca(
  counts_dat,
  sample_metadata,
  sample_id_colname = NULL,
  feature_id_colname = NULL
)
```

## Arguments

- counts_dat:

  data frame of feature counts (e.g. from the counts slot of a
  `multiOmicDataSet`).

- sample_metadata:

  sample metadata as a data frame or tibble.

- sample_id_colname:

  The column from the sample metadata containing the sample names. The
  names in this column must exactly match the names used as the sample
  column names of your input Counts Matrix. (Default: `NULL` - first
  column in the sample metadata will be used.)

- feature_id_colname:

  The column from the counts dataa containing the Feature IDs (Usually
  Gene or Protein ID). This is usually the first column of your input
  Counts Matrix. Only columns of Text type from your input Counts Matrix
  will be available to select for this parameter. (Default: `NULL` -
  first column in the counts matrix will be used.)

## Value

data frame with statistics for each principal component

## See also

Other PCA functions:
[`plot_pca()`](https://ccbr.github.io/MOSuite/reference/plot_pca.md),
[`plot_pca_2d()`](https://ccbr.github.io/MOSuite/reference/plot_pca_2d.md),
[`plot_pca_3d()`](https://ccbr.github.io/MOSuite/reference/plot_pca_3d.md)

## Examples

``` r
calc_pca(nidap_raw_counts, nidap_sample_metadata) %>% head()
#> # A tibble: 6 × 10
#>   Sample    PC  value std.dev percent cumulative Group Replicate Batch Label
#>   <chr>  <dbl>  <dbl>   <dbl>   <dbl>      <dbl> <chr>     <dbl> <dbl> <chr>
#> 1 A1         1 -40.6     61.8   21.2       0.212 A             1     1 A1   
#> 2 A1         2  25.2     56.0   17.4       0.386 A             1     1 A1   
#> 3 A1         3  -9.11    48.5   13.1       0.517 A             1     1 A1   
#> 4 A1         4  21.7     46.9   12.2       0.639 A             1     1 A1   
#> 5 A1         5 -29.3     42.6   10.1       0.740 A             1     1 A1   
#> 6 A1         6 -74.5     41.7    9.68      0.837 A             1     1 A1   
```
