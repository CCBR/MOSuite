# Get vector of colors for observations in one column of a data frame

Get vector of colors for observations in one column of a data frame

## Usage

``` r
get_colors_vctr(dat, colname, palette_fun = select_mosuite_colors, ...)
```

## Arguments

- dat:

  data frame

- colname:

  column name in `dat`

- palette_fun:

  Function for selecting colors. Assumed to contain `n` for the number
  of colors. Defaults to MOSuite's default plot palette. To use the
  previous R default palette behavior, pass
  [`grDevices::palette.colors`](https://rdrr.io/r/grDevices/palette.html).

- ...:

  additional arguments forwarded to `palette_fun`

## Value

named vector of colors for each unique observation in `dat[[colname]]`.
Factor columns use factor level order; other columns use first-observed
order.
