# Get constructor-style default colors for multiple multiOmicDataSet columns

Returns a named list of colors, one entry per requested column. Stored
colors are preferred and missing entries fall back to constructor-style
defaults.

## Usage

``` r
get_moo_default_color_list(moo, colnames, palette = mosuite_palette)
```

## Arguments

- moo:

  `multiOmicDataSet` object

- colnames:

  character vector of column names in `moo@sample_meta`

- palette:

  Character vector of colors to assign. Defaults to `mosuite_palette`.

## Value

Named list of color vectors.
