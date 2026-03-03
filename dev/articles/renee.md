# RSEM counts from RENEE

``` r
library(MOSuite)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

## RENEE dataset

``` r
# replace these lines with the actual paths to your files
gene_counts_tsv <- system.file("extdata",
  "RSEM.genes.expected_count.all_samples.txt.gz",
  package = "MOSuite"
)
metadata_tsv <- system.file("extdata", "sample_metadata.tsv.gz",
  package = "MOSuite"
)

# create multi-omic object
moo <- create_multiOmicDataSet_from_files(
  sample_meta_filepath = metadata_tsv,
  feature_counts_filepath = gene_counts_tsv
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
head(moo@sample_meta)
#> # A tibble: 4 × 2
#>   sample_id condition
#>   <chr>     <chr>    
#> 1 KO_S3     knockout 
#> 2 KO_S4     knockout 
#> 3 WT_S1     wildtype 
#> 4 WT_S2     wildtype
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
```

``` r
moo <- moo |>
  clean_raw_counts() |>
  filter_counts(
    group_colname = "condition",
    label_colname = "sample_id",
    minimum_count_value_to_be_considered_nonzero = 1,
    minimum_number_of_samples_with_nonzero_counts_in_total = 1,
    minimum_number_of_samples_with_nonzero_counts_in_a_group = 1,
  ) |>
  normalize_counts(
    group_colname = "condition",
    label_colname = "sample_id"
  ) |>
  diff_counts(
    covariates_colnames = "condition",
    contrast_colname = "condition",
    contrasts = c("knockout-wildtype")
  ) |>
  filter_diff(
    significance_cutoff = 0.05,
    significance_column = "adjpval",
    change_column = "logFC",
    change_cutoff = 1
  )
#> Saving 7.29 x 4.51 in image
#> * cleaning raw counts
#> 
#> Not able to identify multiple id's in gene_id
#> 
#> Columns that can be used to aggregate gene information gene_id
#> 
#> Aggregating the counts for the same ID in different chromosome locations.
#> Column used to Aggregate duplicate IDs: gene_id
#> Number of rows before Collapse: 58929
#> 
#> no duplicated IDs in gene_id
#> 
#> * filtering clean counts
#> 
#> Number of features after filtering: 291
#> 
#> colors_for_plots NULL
#> 
#> colors_for_plots character
#> 
#> Saving 7.29 x 4.51 in image
#> Saving 7.29 x 4.51 in image
#> * normalizing filt counts
#> 
#> Total number of features included: 291
#> 
#> Saving 7.29 x 4.51 in image
#> Saving 7.29 x 4.51 in image
#> Sample columns: KO_S3, Sample columns: KO_S4, Sample columns: WT_S1, Sample columns: WT_S2
#> 
#> * differential counts
#> 
#> Setting first column of `counts` as gene annotation.
#> 
#> Total number of genes included: 291
#> 
#> Saving 7.29 x 4.51 in image
#> `geom_smooth()` using method = 'loess' and formula = 'y ~ x'
#> * filtering differential features
#> 
#> Total number of genes selected with adjpval < 0.05 and | logFC | ≥ 1 is sum(selgenes)
#> 
#> Saving 7.29 x 4.51 in image

moo@counts$norm$voom |> head()
#>              gene_id     KO_S3     KO_S4     WT_S1     WT_S2
#> 1  ENSG00000215458.8 11.075196 12.348091  8.816153 10.004874
#> 2 ENSG00000160179.18  9.608634 12.770317 12.348091 12.236996
#> 3  ENSG00000258017.1  9.608634  8.816153  8.816153  8.816153
#> 4  ENSG00000282393.1  8.816153  9.608634  8.816153  8.816153
#> 5  ENSG00000286104.1  9.608634  8.816153  8.816153  8.816153
#> 6  ENSG00000274422.1  8.816153  9.608634  8.816153  8.816153
```

## The multiOmicDataSet object structure

``` r
str(moo)
#> <MOSuite::multiOmicDataSet>
#>  @ sample_meta: spc_tbl_ [4 × 2] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
#>  $ sample_id: chr [1:4] "KO_S3" "KO_S4" "WT_S1" "WT_S2"
#>  $ condition: chr [1:4] "knockout" "knockout" "wildtype" "wildtype"
#>  - attr(*, "spec")=
#>   .. cols(
#>   ..   sample_id = col_character(),
#>   ..   condition = col_character()
#>   .. )
#>  - attr(*, "problems")=<externalptr> 
#>  @ annotation : tibble [58,929 × 2] (S3: tbl_df/tbl/data.frame)
#>  $ gene_id : chr [1:58929] "ENSG00000121410.11" "ENSG00000268895.5" "ENSG00000148584.15" "ENSG00000175899.14" ...
#>  $ GeneName: chr [1:58929] "A1BG" "A1BG-AS1" "A1CF" "A2M" ...
#>  @ counts     :List of 4
#>  .. $ raw  : tibble [58,929 × 5] (S3: tbl_df/tbl/data.frame)
#>  ..  ..$ gene_id: chr [1:58929] "ENSG00000121410.11" "ENSG00000268895.5" "ENSG00000148584.15" "ENSG00000175899.14" ...
#>  ..  ..$ KO_S3  : num [1:58929] 0 0 0 0 0 0 0 0 0 0 ...
#>  ..  ..$ KO_S4  : num [1:58929] 0 0 0 0 0 0 0 0 0 0 ...
#>  ..  ..$ WT_S1  : num [1:58929] 0 0 0 0 0 0 0 0 0 0 ...
#>  ..  ..$ WT_S2  : num [1:58929] 0 0 0 0 0 0 0 0 0 0 ...
#>  .. $ clean:'data.frame':    58929 obs. of  5 variables:
#>  ..  ..$ gene_id: chr [1:58929] "ENSG00000121410.11" "ENSG00000268895.5" "ENSG00000148584.15" "ENSG00000175899.14" ...
#>  ..  ..$ KO_S3  : num [1:58929] 0 0 0 0 0 0 0 0 0 0 ...
#>  ..  ..$ KO_S4  : num [1:58929] 0 0 0 0 0 0 0 0 0 0 ...
#>  ..  ..$ WT_S1  : num [1:58929] 0 0 0 0 0 0 0 0 0 0 ...
#>  ..  ..$ WT_S2  : num [1:58929] 0 0 0 0 0 0 0 0 0 0 ...
#>  .. $ filt :'data.frame':    291 obs. of  5 variables:
#>  ..  ..$ gene_id: chr [1:291] "ENSG00000215458.8" "ENSG00000160179.18" "ENSG00000258017.1" "ENSG00000282393.1" ...
#>  ..  ..$ KO_S3  : num [1:291] 2 1 1 0 1 0 0 0 3 0 ...
#>  ..  ..$ KO_S4  : num [1:291] 4 5 0 1 0 1 0 0 3 0 ...
#>  ..  ..$ WT_S1  : num [1:291] 0 6 0 0 0 0 46 33 9 0 ...
#>  ..  ..$ WT_S2  : num [1:291] 1 7 0 0 0 0 80 31 11 1 ...
#>  .. $ norm :List of 1
#>  ..  ..$ voom:'data.frame':  291 obs. of  5 variables:
#>  ..  .. ..$ gene_id: chr [1:291] "ENSG00000215458.8" "ENSG00000160179.18" "ENSG00000258017.1" "ENSG00000282393.1" ...
#>  ..  .. ..$ KO_S3  : num [1:291] 11.08 9.61 9.61 8.82 9.61 ...
#>  ..  .. ..$ KO_S4  : num [1:291] 12.35 12.77 8.82 9.61 8.82 ...
#>  ..  .. ..$ WT_S1  : num [1:291] 8.82 12.35 8.82 8.82 8.82 ...
#>  ..  .. ..$ WT_S2  : num [1:291] 10 12.24 8.82 8.82 8.82 ...
#>  @ analyses   :List of 3
#>  .. $ colors   :List of 2
#>  ..  ..$ sample_id: Named chr [1:4] "#000000" "#E69F00" "#56B4E9" "#009E73"
#>  ..  .. ..- attr(*, "names")= chr [1:4] "KO_S3" "KO_S4" "WT_S1" "WT_S2"
#>  ..  ..$ condition: Named chr [1:2] "#000000" "#E69F00"
#>  ..  .. ..- attr(*, "names")= chr [1:2] "knockout" "wildtype"
#>  .. $ diff     :List of 1
#>  ..  ..$ knockout-wildtype:'data.frame': 291 obs. of  6 variables:
#>  ..  .. ..$ gene_id: chr [1:291] "ENSG00000215458.8" "ENSG00000160179.18" "ENSG00000258017.1" "ENSG00000282393.1" ...
#>  ..  .. ..$ FC     : num [1:291] 5.02 -2.14 1.32 1.32 1.32 ...
#>  ..  .. ..$ logFC  : num [1:291] 2.327 -1.094 0.396 0.396 0.396 ...
#>  ..  .. ..$ tstat  : num [1:291] 3.098 -0.948 0.945 0.945 0.945 ...
#>  ..  .. ..$ pval   : num [1:291] 0.0596 0.4186 0.4196 0.4196 0.4196 ...
#>  ..  .. ..$ adjpval: num [1:291] 0.175 0.49 0.49 0.49 0.49 ...
#>  .. $ diff_filt:'data.frame':    58 obs. of  6 variables:
#>  ..  ..$ gene_id                  : chr [1:58] "ENSG00000154734.15" "ENSG00000154736.6" "ENSG00000232855.6" "ENSG00000231324.1" ...
#>  ..  ..$ knockout-wildtype_FC     : num [1:58] -54.2 -35.6 -7.21 -2.76 6.7 -28.4 -6.43 -27.8 -5.6 6.25 ...
#>  ..  ..$ knockout-wildtype_logFC  : num [1:58] -5.76 -5.15 -2.85 -1.46 2.74 -4.83 -2.69 -4.8 -2.49 2.64 ...
#>  ..  ..$ knockout-wildtype_tstat  : num [1:58] -89.9 -26.5 -13.6 -6.39 10.5 -7.27 -12.4 -15.7 -10.4 8.19 ...
#>  ..  ..$ knockout-wildtype_pval   : num [1:58] 7.21e-06 2.10e-04 1.30e-03 9.94e-03 2.67e-03 7.08e-03 1.68e-03 8.81e-04 2.72e-03 5.14e-03 ...
#>  ..  ..$ knockout-wildtype_adjpval: num [1:58] 0.00119 0.00595 0.0145 0.0402 0.0208 0.0326 0.0158 0.0117 0.0208 0.0245 ...
```
