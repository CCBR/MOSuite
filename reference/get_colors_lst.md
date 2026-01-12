# Create named list of default colors for plotting

Create named list of default colors for plotting

## Usage

``` r
get_colors_lst(sample_metadata, palette_fun = grDevices::palette.colors, ...)
```

## Arguments

- sample_metadata:

  sample metadata as a data frame or tibble. The first column is assumed
  to contain the sample IDs which must correspond to column names in the
  raw counts.

- palette_fun:

  Function for selecting colors. Assumed to contain `n` for the number
  of colors. Default:
  [`grDevices::palette.colors()`](https://rdrr.io/r/grDevices/palette.html)

- ...:

  additional arguments forwarded to `palette_fun`

## Value

named list, with each column in `sample_metadata` containing entry with
a named vector of colors

## Examples

``` r
get_colors_lst(nidap_sample_metadata)
#> $Sample
#>        A1        A2        A3        B1        B2        B3        C1        C2 
#> "#000000" "#E69F00" "#56B4E9" "#009E73" "#F0E442" "#0072B2" "#D55E00" "#CC79A7" 
#>        C3 
#> "#999999" 
#> 
#> $Group
#>         A         B         C 
#> "#000000" "#E69F00" "#56B4E9" 
#> 
#> $Replicate
#>         1         2         3 
#> "#000000" "#E69F00" "#56B4E9" 
#> 
#> $Batch
#>         1         2 
#> "#000000" "#E69F00" 
#> 
#> $Label
#>        A1        A2        A3        B1        B2        B3        C1        C2 
#> "#000000" "#E69F00" "#56B4E9" "#009E73" "#F0E442" "#0072B2" "#D55E00" "#CC79A7" 
#>        C3 
#> "#999999" 
#> 
if (FALSE) { # \dontrun{
get_colors_lst(nidap_sample_metadata, palette_fun = RColorBrewer::brewer.pal, name = "Set3")
} # }
```
