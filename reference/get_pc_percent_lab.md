# Get label for Principal Component with percent of variation

Get label for Principal Component with percent of variation

## Usage

``` r
get_pc_percent_lab(pca_df, pc)
```

## Arguments

- pca_df:

  data frame from
  [`calc_pca()`](https://ccbr.github.io/MOSuite/reference/calc_pca.md)

- pc:

  which principal component to report (e.g. `1`)

## Value

glue string formatted with PC's percent of variation

## Examples

``` r
if (FALSE) { # \dontrun{
data.frame(PC = c(1, 2, 3), percent = c(40, 10, 0.5)) |>
  get_pc_percent_lab(2)
} # }
```
