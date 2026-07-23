# Plot read depth for `data.frame`

Plot read depth for `data.frame`

## Arguments

- sample_metadata:

  sample metadata dataframe, required when `group_colname` is supplied.

- sample_id_colname:

  column in sample metadata containing sample IDs.

- group_colname:

  sample metadata column used to color bars. Leave blank to use the
  current single-color bar fill.

- color_values:

  colors used when `group_colname` is supplied. Named vectors are
  matched to group values; unnamed vectors follow group order and are
  extended with MOSuite colors when too few colors are supplied.
  Defaults to `NULL`; when `NULL`, `mosuite_palette` is used.

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
