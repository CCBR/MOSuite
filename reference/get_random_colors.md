# Get random colors.

Note: this function is not guaranteed to create a color blind friendly
palette. Consider using other palettes such as
`RColorBrewer::display.brewer.all(colorblindFriendly = TRUE)`.

## Usage

``` r
get_random_colors(num_colors, n = 2000)
```

## Arguments

- num_colors:

  number of colors to select.

- n:

  number of random RGB values to generate in the color space.

## Value

vector of random colors in hex format.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(10)
get_random_colors(5)
} # }
```
