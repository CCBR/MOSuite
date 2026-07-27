# Default MOSuite color palette

A character vector of 12 hex color codes used as the default palette
throughout MOSuite plots. Colors are assigned to groups in the order
they appear (or by factor level order). Use
[`set_color_pal()`](https://ccbr.github.io/MOSuite/reference/set_color_pal.md)
to override the palette for a specific metadata column.

## Usage

``` r
mosuite_palette
```

## Format

A character vector of length 12.

## See also

[`select_mosuite_colors()`](https://ccbr.github.io/MOSuite/reference/select_mosuite_colors.md),
[`set_color_pal()`](https://ccbr.github.io/MOSuite/reference/set_color_pal.md)

## Examples

``` r
mosuite_palette
#>  [1] "#5954d6" "#e1562c" "#b80058" "#00c6f8" "#d163e6" "#00a76c" "#ff9287"
#>  [8] "#008cf9" "#006e00" "#796880" "#FFA500" "#878500"
scales::show_col(mosuite_palette)
```
