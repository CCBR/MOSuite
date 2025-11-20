# Plot read depth for multiOmicDataSet

Plot read depth for multiOmicDataSet

## Arguments

- moo_counts:

  `multiOmicDataSet` containing `count_type` & `sub_count_type` in the
  counts slot

- count_type:

  the type of counts to use. Must be a name in the counts slot
  (`names(moo@counts)`).

- sub_count_type:

  used if `count_type` is a list in the counts slot: specify the sub
  count type within the list. Must be a name in
  `names(moo@counts[[count_type]])`.

## Value

ggplot barplot

## See also

[plot_read_depth](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.md)
generic

Other plotters for multiOmicDataSets:
[`plot_corr_heatmap_moo`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap_moo.md),
[`plot_expr_heatmap_moo`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap_moo.md),
[`plot_histogram_moo`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram_moo.md),
[`plot_pca_moo`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_moo.md)

## Examples

``` r
# multiOmicDataSet
moo <- multiOmicDataSet(
  sample_metadata = nidap_sample_metadata,
  anno_dat = data.frame(),
  counts_lst = list(
    "raw" = nidap_raw_counts,
    "clean" = nidap_clean_raw_counts
  )
)

plot_read_depth(moo, count_type = "clean")

```
