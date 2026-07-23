# Resolve plotting colors for one column

Uses `color_values` when supplied; otherwise generates colors with
[`get_colors_vctr()`](https://ccbr.github.io/MOSuite/dev/reference/get_colors_vctr.md).
If too few colors are provided, missing colors are generated and
appended.

## Usage

``` r
resolve_plot_colors(
  dat,
  colname,
  color_values = NULL,
  palette_fun = select_mosuite_colors,
  ...
)
```

## Arguments

- dat:

  data frame

- colname:

  column name in `dat`

- color_values:

  optional vector of colors

- palette_fun:

  function used to generate colors

- ...:

  additional arguments forwarded to `palette_fun`

## Value

named vector of colors matching observed values in `colname`
