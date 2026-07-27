# Display the mosuite color palette

Plots each color in `mosuite_palette` as a labeled tile with its hex
code displayed below. The plot is rendered at a width proportional to
the number of colors so labels remain horizontal and legible.

## Usage

``` r
display_palette(palette = mosuite_palette)
```

## Arguments

- palette:

  Character vector of hex color codes. Defaults to `mosuite_palette`.

## Value

Invisibly returns the underlying
[ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
display_palette()

display_palette(c("#FF0000", "#00FF00", "#0000FF"))
```
