# Create named list of default colors for plotting

Create named list of default colors for plotting

## Usage

``` r
get_colors_lst(sample_metadata, palette = mosuite_palette)
```

## Arguments

- sample_metadata:

  sample metadata as a data frame or tibble. The first column is assumed
  to contain the sample IDs which must correspond to column names in the
  raw counts.

- palette:

  Character vector of colors to assign. Defaults to `mosuite_palette`.

## Value

named list, with each column in `sample_metadata` containing a
corresponding entry with a named vector of colors

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
#> "#00c6f8" "#d163e6" "#00a76c" 
#> 
#> $Batch
#>         1         2 
#> "#ff9287" "#008cf9" 
#> 
#> $Label
#>        A1        A2        A3        B1        B2        B3        C1        C2 
#> "#5954d6" "#e1562c" "#b80058" "#00c6f8" "#d163e6" "#00a76c" "#ff9287" "#008cf9" 
#>        C3 
#> "#006e00" 
#> 
get_colors_lst(nidap_sample_metadata, palette = RColorBrewer::brewer.pal(12, "Set3"))
#> $Sample
#>        A1        A2        A3        B1        B2        B3        C1        C2 
#> "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462" "#B3DE69" "#FCCDE5" 
#>        C3 
#> "#D9D9D9" 
#> 
#> $Group
#>         A         B         C 
#> "#8DD3C7" "#FFFFB3" "#BEBADA" 
#> 
#> $Replicate
#>         1         2         3 
#> "#FB8072" "#80B1D3" "#FDB462" 
#> 
#> $Batch
#>         1         2 
#> "#B3DE69" "#FCCDE5" 
#> 
#> $Label
#>        A1        A2        A3        B1        B2        B3        C1        C2 
#> "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462" "#B3DE69" "#FCCDE5" 
#>        C3 
#> "#D9D9D9" 
#> 
```
