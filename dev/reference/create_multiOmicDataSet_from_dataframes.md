# Construct a multiOmicDataSet object from data frames

Construct a multiOmicDataSet object from data frames

## Usage

``` r
create_multiOmicDataSet_from_dataframes(
  sample_metadata,
  counts_dat,
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  count_type = "raw"
)
```

## Arguments

- sample_metadata:

  sample metadata as a data frame or tibble. The first column is assumed
  to contain the sample IDs which must correspond to column names in the
  raw counts.

- counts_dat:

  data frame of feature counts (e.g. expected feature counts from RSEM).

- sample_id_colname:

  name of the column in `sample_metadata` that contains the sample IDs.
  (Default: `NULL` - first column in the sample metadata will be used.)

- feature_id_colname:

  name of the column in `counts_dat` that contains feature/gene IDs.
  (Default: `NULL` - first column in the count data will be used.)

- count_type:

  type to assign the values of `counts_dat` to in the `counts` slot

## Value

[multiOmicDataSet](https://ccbr.github.io/MOSuite/dev/reference/multiOmicDataSet.md)
object

## See also

Other moo constructors:
[`create_multiOmicDataSet_from_files()`](https://ccbr.github.io/MOSuite/dev/reference/create_multiOmicDataSet_from_files.md),
[`multiOmicDataSet()`](https://ccbr.github.io/MOSuite/dev/reference/multiOmicDataSet.md)

## Examples

``` r
sample_meta <- data.frame(
  sample_id = c("KO_S3", "KO_S4", "WT_S1", "WT_S2"),
  condition = factor(
    c("knockout", "knockout", "wildtype", "wildtype"),
    levels = c("wildtype", "knockout")
  )
)
moo <- create_multiOmicDataSet_from_dataframes(sample_meta, gene_counts)
head(moo@sample_meta)
#>   sample_id condition
#> 1     KO_S3  knockout
#> 2     KO_S4  knockout
#> 3     WT_S1  wildtype
#> 4     WT_S2  wildtype
head(moo@counts$raw)
#> # A tibble: 6 × 5
#>   gene_id            KO_S3 KO_S4 WT_S1 WT_S2
#>   <chr>              <dbl> <dbl> <dbl> <dbl>
#> 1 ENSG00000121410.11     0     0     0     0
#> 2 ENSG00000268895.5      0     0     0     0
#> 3 ENSG00000148584.15     0     0     0     0
#> 4 ENSG00000175899.14     0     0     0     0
#> 5 ENSG00000245105.3      0     0     0     0
#> 6 ENSG00000166535.20     0     0     0     0
head(moo@annotation)
#> # A tibble: 6 × 2
#>   gene_id            GeneName
#>   <chr>              <chr>   
#> 1 ENSG00000121410.11 A1BG    
#> 2 ENSG00000268895.5  A1BG-AS1
#> 3 ENSG00000148584.15 A1CF    
#> 4 ENSG00000175899.14 A2M     
#> 5 ENSG00000245105.3  A2M-AS1 
#> 6 ENSG00000166535.20 A2ML1   

sample_meta_nidap <- readr::read_csv(system.file("extdata", "nidap",
  "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz",
  package = "MOSuite"
))
#> Rows: 9 Columns: 5
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (3): Sample, Group, Label
#> dbl (2): Replicate, Batch
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
raw_counts_nidap <- readr::read_csv(system.file("extdata", "nidap", "Raw_Counts.csv.gz",
  package = "MOSuite"
))
#> Rows: 43280 Columns: 10
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (1): GeneName
#> dbl (9): A1, A2, A3, B1, B2, B3, C1, C2, C3
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
moo_nidap <- create_multiOmicDataSet_from_dataframes(sample_meta_nidap, raw_counts_nidap)
```
