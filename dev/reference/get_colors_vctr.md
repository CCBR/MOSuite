# Get vector of colors for observations in one column of a data frame

Assigns one color per unique observed value in `dat[[colname]]`, drawn
from `palette` starting at `color_offset`. If the palette is too short,
falls back to
[`get_random_colors()`](https://ccbr.github.io/MOSuite/dev/reference/get_random_colors.md).
Factor columns use factor-level order; other columns use first-observed
order.

## Usage

``` r
get_colors_vctr(dat, colname, palette = mosuite_palette, color_offset = 0L)
```

## Arguments

- dat:

  data frame

- colname:

  column name in `dat`

- palette:

  Character vector of colors to assign. Defaults to `mosuite_palette`.

- color_offset:

  integer; number of palette colors to skip before assigning colors to
  this column's values. Used by
  [`get_colors_lst()`](https://ccbr.github.io/MOSuite/dev/reference/get_colors_lst.md)
  to avoid repeating colors across columns with few unique values.

## Value

Named character vector of hex colors, one per unique observed value in
`dat[[colname]]`.

## Examples

``` r
get_colors_vctr(nidap_sample_metadata, "Group")
#>         A         B         C 
#> "#5954d6" "#e1562c" "#b80058" 
get_colors_vctr(nidap_sample_metadata, "Group", color_offset = 3L)
#>         A         B         C 
#> "#00c6f8" "#d163e6" "#00a76c" 
```
