# Plot 2D or 3D PCA for multiOmicDataset

Plot 2D or 3D PCA for multiOmicDataset

## Arguments

- count_type:

  the type of counts to use. Must be a name in the counts slot
  (`names(moo@counts)`).

- sub_count_type:

  used if `count_type` is a list in the counts slot: specify the sub
  count type within the list. Must be a name in
  `names(moo@counts[[count_type]])`.

## Value

PCA plot

## See also

[`plot_pca()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md)
generic

Other plotters for multiOmicDataSets:
[`plot_corr_heatmap,MOSuite::multiOmicDataSet-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap-multiOmicDataSet.md),
[`plot_histogram,MOSuite::multiOmicDataSet-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram.multiOmicDataSet.md),
[`plot_read_depth,MOSuite::multiOmicDataSet-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.multiOmicDataSet.md)
