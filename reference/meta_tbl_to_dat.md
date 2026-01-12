# Convert sample metadata from a tibble to a dataframe with sample IDs as row names

Convert sample metadata from a tibble to a dataframe with sample IDs as
row names

## Usage

``` r
meta_tbl_to_dat(meta_tbl, sample_id_colname = sample_id)
```

## Arguments

- meta_tbl:

  tibble with `sample_id` column

- sample_id_colname:

  name of the column in `sample_metadata` that contains the sample IDs.
  (Default: `NULL` - first column in the sample metadata will be used.)

## Value

dataframe where row names are the sample IDs

## Examples

``` r
if (FALSE) { # \dontrun{
sample_meta_tbl <- readr::read_tsv(system.file("extdata",
  "sample_metadata.tsv.gz",
  package = "MOSuite"
))
head(sample_meta_tbl)
meta_tbl_to_dat(sample_meta_tbl)
} # }
```
