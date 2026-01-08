# Plot correlation heatmap for multiOmicDataSet

Plot correlation heatmap for multiOmicDataSet

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

- ...:

  arguments forwarded to method
  [plot_corr_heatmap_dat](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap_dat.md)

## See also

[plot_corr_heatmap](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap.md)
generic

Other plotters for multiOmicDataSets:
[`plot_histogram_moo`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram_moo.md),
[`plot_pca_moo`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_moo.md),
[`plot_read_depth_moo`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth_moo.md)
