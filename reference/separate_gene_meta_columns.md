# Separate gene metadata column

Separate gene metadata column

## Usage

``` r
separate_gene_meta_columns(counts_dat, split_gene_name = TRUE)
```

## Arguments

- counts_dat:

  dataframe with raw counts data

- split_gene_name:

  If `TRUE`, split the gene name column by any of these special
  characters: `,|_-:`

## Value

dataframe with metadata separated
