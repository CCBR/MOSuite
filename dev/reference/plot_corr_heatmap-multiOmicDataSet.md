# Plot correlation heatmap for multiOmicDataSet

Plot correlation heatmap for multiOmicDataSet

## Arguments

- moo_counts:

  a `multiOmicDataSet` object

- count_type:

  the type of counts to use. Must be a name in the counts slot
  (`names(moo@counts)`).

- sub_count_type:

  used if `count_type` is a list in the counts slot: specify the sub
  count type within the list. Must be a name in
  `names(moo@counts[[count_type]])`.

- ...:

  additional arguments forwarded to
  [`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap.md)
  for `data.frame`

## See also

[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap.md)
generic

Other plotters for multiOmicDataSets:
[`plot_histogram,MOSuite::multiOmicDataSet-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram.multiOmicDataSet.md),
[`plot_pca,MOSuite::multiOmicDataSet-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.multiOmicDataSet.md),
[`plot_read_depth,MOSuite::multiOmicDataSet-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.multiOmicDataSet.md)
