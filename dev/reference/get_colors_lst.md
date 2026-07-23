# Create named list of default colors for plotting

Create named list of default colors for plotting

## Usage

``` r
get_colors_lst(sample_metadata, palette_fun = select_mosuite_colors, ...)
```

## Arguments

- sample_metadata:

  sample metadata as a data frame or tibble. The first column is assumed
  to contain the sample IDs which must correspond to column names in the
  raw counts.

- palette_fun:

  Function for selecting colors. Assumed to contain `n` for the number
  of colors. Defaults to MOSuite's default plot palette. To use the
  previous R default palette behavior, pass
  [`grDevices::palette.colors`](https://rdrr.io/r/grDevices/palette.html).

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
#> "#5954d6" "#e1562c" "#b80058" "#00c6f8" "#d163e6" "#00a76c" "#ff9287" "#008cf9" 
#>        C3 
#> "#006e00" 
#> 
#> $Group
#>         A         B         C 
#> "#5954d6" "#e1562c" "#b80058" 
#> 
#> $Replicate
#>         1         2         3 
#> "#5954d6" "#e1562c" "#b80058" 
#> 
#> $Batch
#>         1         2 
#> "#5954d6" "#e1562c" 
#> 
#> $Label
#>        A1        A2        A3        B1        B2        B3        C1        C2 
#> "#5954d6" "#e1562c" "#b80058" "#00c6f8" "#d163e6" "#00a76c" "#ff9287" "#008cf9" 
#>        C3 
#> "#006e00" 
#> 
if (FALSE) { # \dontrun{
get_colors_lst(nidap_sample_metadata, palette_fun = RColorBrewer::brewer.pal, name = "Set3")
} # }
```
