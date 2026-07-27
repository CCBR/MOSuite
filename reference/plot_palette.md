# Plot a palette tile strip

Renders a data frame with columns `hex` and `idx` as a row of colored
tiles, each labeled with its hex code. Used internally by
[`display_palette()`](https://ccbr.github.io/MOSuite/reference/display_palette.md)
and
[`display_colors()`](https://ccbr.github.io/MOSuite/reference/display_colors.md).

## Usage

``` r
plot_palette(dat)
```

## Arguments

- dat:

  data frame with columns `hex` (hex color codes) and `idx` (factor,
  used for faceting)

## Value

a [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object
