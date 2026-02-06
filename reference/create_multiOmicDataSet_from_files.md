# Construct a multiOmicDataSet object from text files (e.g. TSV, CSV).

Construct a multiOmicDataSet object from text files (e.g. TSV, CSV).

## Usage

``` r
create_multiOmicDataSet_from_files(
  sample_meta_filepath,
  feature_counts_filepath,
  count_type = "raw",
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  delim = NULL,
  ...
)
```

## Arguments

- sample_meta_filepath:

  path to text file with sample IDs and metadata for differential
  analysis.

- feature_counts_filepath:

  path to text file of expected feature counts (e.g. gene counts from
  RSEM).

- count_type:

  type to assign the values of `counts_dat` to in the `counts` slot

- sample_id_colname:

  name of the column in `sample_metadata` that contains the sample IDs.
  (Default: `NULL` - first column in the sample metadata will be used.)

- feature_id_colname:

  name of the column in `counts_dat` that contains feature/gene IDs.
  (Default: `NULL` - first column in the count data will be used.)

- delim:

  Delimiter used in the input files. Any delimiter accepted by
  [`readr::read_delim()`](https://readr.tidyverse.org/reference/read_delim.html)
  can be used. If the files are in CSV format, set `delim = ','`; for
  TSV format, set `delim = '\t'`.

- ...:

  additional arguments forwarded to
  [`readr::read_delim()`](https://readr.tidyverse.org/reference/read_delim.html).

## Value

[multiOmicDataSet](https://ccbr.github.io/MOSuite/reference/multiOmicDataSet.md)
object

## See also

Other moo constructors:
[`create_multiOmicDataSet_from_dataframes()`](https://ccbr.github.io/MOSuite/reference/create_multiOmicDataSet_from_dataframes.md),
[`multiOmicDataSet()`](https://ccbr.github.io/MOSuite/reference/multiOmicDataSet.md)

## Examples

``` r
moo <- create_multiOmicDataSet_from_files(
  sample_meta_filepath = system.file("extdata",
    "sample_metadata.tsv.gz",
    package = "MOSuite"
  ),
  feature_counts_filepath = system.file("extdata",
    "RSEM.genes.expected_count.all_samples.txt.gz",
    package = "MOSuite"
  ),
  delim = "\t"
)
#> Rows: 58929 Columns: 6
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (2): gene_id, GeneName
#> dbl (4): KO_S3, KO_S4, WT_S1, WT_S2
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 4 Columns: 2
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (2): sample_id, condition
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
moo@counts$raw |> head()
#> # A tibble: 6 × 5
#>   gene_id            KO_S3 KO_S4 WT_S1 WT_S2
#>   <chr>              <dbl> <dbl> <dbl> <dbl>
#> 1 ENSG00000121410.11     0     0     0     0
#> 2 ENSG00000268895.5      0     0     0     0
#> 3 ENSG00000148584.15     0     0     0     0
#> 4 ENSG00000175899.14     0     0     0     0
#> 5 ENSG00000245105.3      0     0     0     0
#> 6 ENSG00000166535.20     0     0     0     0
moo@sample_meta
#> # A tibble: 4 × 2
#>   sample_id condition
#>   <chr>     <chr>    
#> 1 KO_S3     knockout 
#> 2 KO_S4     knockout 
#> 3 WT_S1     wildtype 
#> 4 WT_S2     wildtype 

moo_nidap <- create_multiOmicDataSet_from_files(
  system.file("extdata", "nidap",
    "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz",
    package = "MOSuite"
  ),
  system.file("extdata", "nidap", "Raw_Counts.csv.gz", package = "MOSuite"),
  delim = ","
)
#> Rows: 43280 Columns: 10
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (1): GeneName
#> dbl (9): A1, A2, A3, B1, B2, B3, C1, C2, C3
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 9 Columns: 5
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (3): Sample, Group, Label
#> dbl (2): Replicate, Batch
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```
