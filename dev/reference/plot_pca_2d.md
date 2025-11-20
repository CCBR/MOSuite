# Perform and plot a 2D Principal Components Analysis

Perform and plot a 2D Principal Components Analysis

## Usage

``` r
plot_pca_2d(
  counts_dat,
  sample_metadata,
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  group_colname = "Group",
  label_colname = "Label",
  samples_to_rename = NULL,
  color_values = c("#5954d6", "#e1562c", "#b80058", "#00c6f8", "#d163e6", "#00a76c",
    "#ff9287", "#008cf9", "#006e00", "#796880", "#FFA500", "#878500"),
  principal_components = c(1, 2),
  legend_position = "top",
  point_size = 1,
  add_label = TRUE,
  label_font_size = 3,
  label_offset_x_ = 2,
  label_offset_y_ = 2,
  interactive_plots = FALSE
)
```

## Arguments

- counts_dat:

  data frame of feature counts (e.g. expected feature counts from RSEM).

- sample_metadata:

  sample metadata as a data frame or tibble.

- sample_id_colname:

  The column from the sample metadata containing the sample names. The
  names in this column must exactly match the names used as the sample
  column names of your input Counts Matrix. (Default: `NULL` - first
  column in the sample metadata will be used.)

- feature_id_colname:

  The column from the counts dataa containing the Feature IDs (Usually
  Gene or Protein ID). This is usually the first column of your input
  Counts Matrix. Only columns of Text type from your input Counts Matrix
  will be available to select for this parameter. (Default: `NULL` -
  first column in the counts matrix will be used.)

- group_colname:

  The column from the sample metadata containing the sample group
  information. This is usually a column showing to which experimental
  treatments each sample belongs (e.g. WildType, Knockout, Tumor,
  Normal, Before, After, etc.).

- label_colname:

  The column from the sample metadata containing the sample labels as
  you wish them to appear in the plots produced by this template. This
  can be the same Sample Names Column. However, you may desire different
  labels to display on your figure (e.g. shorter labels are sometimes
  preferred on plots). In that case, select the column with your
  preferred Labels here. The selected column should contain unique names
  for each sample. (Default: `NULL` – `sample_id_colname` will be used.)

- samples_to_rename:

  If you do not have a Plot Labels Column in your sample metadata table,
  you can use this parameter to rename samples manually for display on
  the PCA plot. Use "Add item" to add each additional sample for
  renaming. Use the following format to describe which old name (in your
  sample metadata table) you want to rename to which new name: old_name:
  new_name

- color_values:

  vector of colors as hex values or names recognized by R

- principal_components:

  vector with numbered principal components to plot (Default: `c(1,2)`)

- legend_position:

  passed to in `legend.position`
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)

- point_size:

  size for
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)

- add_label:

  whether to add text labels for the points

- label_font_size:

  label font size for the PCA plot

- label_offset_x\_:

  label offset x for the PCA plot

- label_offset_y\_:

  label offset y for the PCA plot

- interactive_plots:

  set to TRUE to make PCA and Histogram plots interactive with `plotly`,
  allowing you to hover your mouse over a point or line to view sample
  information. The similarity heat map will not display if this toggle
  is set to TRUE. Default is FALSE.

## Value

ggplot object

## See also

[plot_pca](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md)
generic

Other PCA functions:
[`calc_pca()`](https://ccbr.github.io/MOSuite/dev/reference/calc_pca.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md),
[`plot_pca_3d()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_3d.md)
