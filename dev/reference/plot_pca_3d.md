# 3D PCA for counts dataframe

3D PCA for counts dataframe

## Usage

``` r
plot_pca_3d(
  counts_dat,
  sample_metadata,
  sample_id_colname = NULL,
  samples_to_rename = NULL,
  group_colname = "Group",
  label_colname = "Label",
  principal_components = c(1, 2, 3),
  point_size = 8,
  label_font_size = 24,
  color_values = c("#5954d6", "#e1562c", "#b80058", "#00c6f8", "#d163e6", "#00a76c",
    "#ff9287", "#008cf9", "#006e00", "#796880", "#FFA500", "#878500"),
  plot_title = "PCA 3D"
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

- samples_to_rename:

  If you do not have a Plot Labels Column in your sample metadata table,
  you can use this parameter to rename samples manually for display on
  the PCA plot. Use "Add item" to add each additional sample for
  renaming. Use the following format to describe which old name (in your
  sample metadata table) you want to rename to which new name: old_name:
  new_name

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

- principal_components:

  vector with numbered principal components to plot (Default:
  `c(1,2,3)`)

- point_size:

  size for
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)

- label_font_size:

  label font size for the PCA plot

- color_values:

  vector of colors as hex values or names recognized by R

- plot_title:

  title for the plot

## Value

[`plotly::plot_ly`](https://rdrr.io/pkg/plotly/man/plot_ly.html) figure

## See also

[plot_pca](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md)
generic

Other PCA functions:
[`calc_pca()`](https://ccbr.github.io/MOSuite/dev/reference/calc_pca.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md),
[`plot_pca_2d()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_2d.md)
