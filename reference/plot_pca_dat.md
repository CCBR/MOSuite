# Plot 2D or 3D PCA for counts dataframe

Plot 2D or 3D PCA for counts dataframe

## Arguments

- moo_counts:

  counts dataframe

- sample_metadata:

  **Required** if `moo_counts` is a `data.frame`: sample metadata as a
  data frame or tibble.

- principal_components:

  vector with numbered principal components to plot. Use 2 for a 2D pca
  with ggplot, or 3 for a 3D pca with plotly. (Default: `c(1,2)`)

- ...:

  additional arguments forwarded to
  [`plot_pca_2d()`](https://ccbr.github.io/MOSuite/reference/plot_pca_2d.md)
  (if 2 PCs) or
  [`plot_pca_3d()`](https://ccbr.github.io/MOSuite/reference/plot_pca_3d.md)
  (if 3 PCs).

## See also

[plot_pca](https://ccbr.github.io/MOSuite/reference/plot_pca.md) generic

Other plotters for counts dataframes:
[`plot_corr_heatmap_dat`](https://ccbr.github.io/MOSuite/reference/plot_corr_heatmap_dat.md),
[`plot_histogram_dat`](https://ccbr.github.io/MOSuite/reference/plot_histogram_dat.md),
[`plot_read_depth_dat`](https://ccbr.github.io/MOSuite/reference/plot_read_depth_dat.md)
