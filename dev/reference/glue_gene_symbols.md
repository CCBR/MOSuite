# Glue gene_id and GeneName columns into one column

Glue gene_id and GeneName columns into one column

## Usage

``` r
glue_gene_symbols(counts_dat)
```

## Arguments

- counts_dat:

  data frame containing gene_id and GeneName columns

## Value

counts_dat with gene_id and GeneName joined with `|` as the new gene_id
column

## Examples

``` r
if (FALSE) { # \dontrun{
gene_counts |>
  glue_gene_symbols() |>
  head()
} # }
```
