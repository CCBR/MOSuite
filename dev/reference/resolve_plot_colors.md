# Resolve plotting colors for one column

Uses `color_values` when supplied; otherwise generates colors with
[`get_colors_vctr()`](https://ccbr.github.io/MOSuite/dev/reference/get_colors_vctr.md).
If `color_values` is named and covers all observed values, it is
returned as-is. If too few colors are provided, missing colors are
generated and appended.

## Usage

``` r
resolve_plot_colors(
  dat,
  colname,
  color_values = NULL,
  palette = mosuite_palette
)
```

## Arguments

- dat:

  data frame

- colname:

  column name in `dat`

- color_values:

  optional named or unnamed character vector of colors

- palette:

  character vector of colors used to generate defaults

## Value

named character vector of colors matching observed values in
`dat[[colname]]`
