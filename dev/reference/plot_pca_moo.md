# Plot 2D or 3D PCA for multiOmicDataset

Plot 2D or 3D PCA for multiOmicDataset

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

- principal_components:

  vector with numbered principal components to plot. Use 2 for a 2D pca
  with ggplot, or 3 for a 3D pca with plotly. (Default: `c(1,2)`)

- ...:

  additional arguments forwarded to
  [`plot_pca_2d()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_2d.md)
  (if 2 PCs) or
  [`plot_pca_3d()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_3d.md)
  (if 3 PCs).

## Value

PCA plot

## See also

[plot_pca](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md)
generic

Other plotters for multiOmicDataSets:
[`plot_corr_heatmap_moo`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap_moo.md),
[`plot_expr_heatmap_moo`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap_moo.md),
[`plot_histogram_moo`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram_moo.md),
[`plot_read_depth_moo`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth_moo.md)
