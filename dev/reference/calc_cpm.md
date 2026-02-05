# Calculate counts-per-million (CPM) on raw counts in a multiOmicDataSet

Calculate counts-per-million (CPM) on raw counts in a multiOmicDataSet

## Usage

``` r
calc_cpm(moo, ...)
```

## Arguments

- moo:

  multiOmicDataSet object

- ...:

  additional arguments to pass to edgeR::cpm()

## Value

multiOmicDataSet with cpm-transformed counts

## Examples

``` r
sample_meta <- data.frame(
  sample_id = c("KO_S3", "KO_S4", "WT_S1", "WT_S2"),
  condition = factor(
    c("knockout", "knockout", "wildtype", "wildtype"),
    levels = c("wildtype", "knockout")
  )
)
moo <- create_multiOmicDataSet_from_dataframes(sample_meta, gene_counts) |>
  calc_cpm()
head(moo@counts$cpm)
#>              gene_id KO_S3 KO_S4 WT_S1 WT_S2
#> 1 ENSG00000121410.11     0     0     0     0
#> 2  ENSG00000268895.5     0     0     0     0
#> 3 ENSG00000148584.15     0     0     0     0
#> 4 ENSG00000175899.14     0     0     0     0
#> 5  ENSG00000245105.3     0     0     0     0
#> 6 ENSG00000166535.20     0     0     0     0
```
