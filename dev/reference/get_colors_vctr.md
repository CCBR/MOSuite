# Get vector of colors for observations in one column of a data frame

Get vector of colors for observations in one column of a data frame

## Usage

``` r
get_colors_vctr(dat, colname, palette_fun = grDevices::palette.colors, ...)
```

## Arguments

- dat:

  data frame

- colname:

  column name in `dat`

- palette_fun:

  Function for selecting colors. Assumed to contain `n` for the number
  of colors. Default:
  [`grDevices::palette.colors()`](https://rdrr.io/r/grDevices/palette.html)

- ...:

  additional arguments forwarded to `palette_fun`

## Value

named vector of colors for each unique observation in `dat$colname`
