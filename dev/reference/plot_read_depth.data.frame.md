# Plot read depth for `data.frame`

Plot read depth for `data.frame`

## Arguments

- ...:

  additional arguments (ignored; accepted for compatibility with the moo
  dispatch)

## Value

ggplot barplot

## See also

[`plot_read_depth()`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.md)
generic

Other plotters for counts dataframes:
[`plot_corr_heatmap,data.frame-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap-data.frame.md),
[`plot_histogram,data.frame-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram.data.frame.md),
[`plot_pca,data.frame-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.data.frame.md)

## Examples

``` r
# dataframe
plot_read_depth(nidap_clean_raw_counts)

```
