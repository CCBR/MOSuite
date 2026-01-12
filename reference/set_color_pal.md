# Set color palette for a single group/column

This allows you to set custom palettes individually for groups in the
dataset

## Usage

``` r
set_color_pal(moo, colname, palette_fun = grDevices::palette.colors, ...)
```

## Arguments

- moo:

  `multiOmicDataSet` object (see
  [`create_multiOmicDataSet_from_dataframes()`](https://ccbr.github.io/MOSuite/reference/create_multiOmicDataSet_from_dataframes.md))

- colname:

  group column name to set the palette for

- palette_fun:

  Function for selecting colors. Assumed to contain `n` for the number
  of colors. Default:
  [`grDevices::palette.colors()`](https://rdrr.io/r/grDevices/palette.html)

- ...:

  additional arguments forwarded to `palette_fun`

## Value

`moo` with colors updated at `moo@analyses$colors$colname`

## See also

Other moo methods:
[`batch_correct_counts()`](https://ccbr.github.io/MOSuite/reference/batch_correct_counts.md),
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/reference/clean_raw_counts.md),
[`diff_counts()`](https://ccbr.github.io/MOSuite/reference/diff_counts.md),
[`filter_counts()`](https://ccbr.github.io/MOSuite/reference/filter_counts.md),
[`filter_diff()`](https://ccbr.github.io/MOSuite/reference/filter_diff.md),
[`normalize_counts()`](https://ccbr.github.io/MOSuite/reference/normalize_counts.md),
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_corr_heatmap.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_expr_heatmap.md),
[`plot_histogram()`](https://ccbr.github.io/MOSuite/reference/plot_histogram.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/reference/plot_pca.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/reference/plot_read_depth.md),
[`run_deseq2()`](https://ccbr.github.io/MOSuite/reference/run_deseq2.md)

## Examples

``` r
moo <- create_multiOmicDataSet_from_dataframes(
  sample_metadata = as.data.frame(nidap_sample_metadata),
  counts_dat = as.data.frame(nidap_raw_counts)
)
moo@analyses$colors$Group
#>         A         B         C 
#> "#000000" "#E69F00" "#56B4E9" 
moo %<>% set_color_pal("Group", palette_fun = RColorBrewer::brewer.pal, name = "Set2")
moo@analyses$colors$Group
#>         A         B         C 
#> "#66C2A5" "#FC8D62" "#8DA0CB" 
```
