# Perform and plot a 2D Principal Components Analysis

Perform and plot a 2D Principal Components Analysis

Perform and plot a 2D Principal Components Analysis

## Usage

``` r
plot_pca_2d(
  moo_counts,
  count_type = NULL,
  sub_count_type = NULL,
  sample_metadata = NULL,
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  group_colname = "Group",
  label_colname = "Label",
  samples_to_rename = NULL,
  color_values = NULL,
  principal_components = c(1, 2),
  legend_position = "top",
  point_size = 5,
  legend_font_size = NULL,
  label_font_size = 3,
  label_offset_x_ = 2,
  label_offset_y_ = 2,
  log_transform = FALSE,
  log_transform_pseudocount = 0.5,
  log_transform_base = "ln",
  interactive_plots = FALSE,
  plots_subdir = "pca",
  plot_filename = "pca_2D.png",
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  ...
)

## S7 method for class <MOSuite::multiOmicDataSet>
plot_pca_2d(
  moo_counts,
  count_type = NULL,
  sub_count_type = NULL,
  sample_metadata = NULL,
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  group_colname = "Group",
  label_colname = "Label",
  samples_to_rename = NULL,
  color_values = NULL,
  principal_components = c(1, 2),
  legend_position = "top",
  point_size = 5,
  legend_font_size = NULL,
  label_font_size = 3,
  label_offset_x_ = 2,
  label_offset_y_ = 2,
  log_transform = FALSE,
  log_transform_pseudocount = 0.5,
  log_transform_base = "ln",
  interactive_plots = FALSE,
  plots_subdir = "pca",
  plot_filename = "pca_2D.png",
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  ...
)

## S7 method for class <data.frame>
plot_pca_2d(
  moo_counts,
  count_type = NULL,
  sub_count_type = NULL,
  sample_metadata = NULL,
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  group_colname = "Group",
  label_colname = "Label",
  samples_to_rename = NULL,
  color_values = NULL,
  principal_components = c(1, 2),
  legend_position = "top",
  point_size = 5,
  legend_font_size = NULL,
  label_font_size = 3,
  label_offset_x_ = 2,
  label_offset_y_ = 2,
  log_transform = FALSE,
  log_transform_pseudocount = 0.5,
  log_transform_base = "ln",
  interactive_plots = FALSE,
  plots_subdir = "pca",
  plot_filename = "pca_2D.png",
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  ...
)
```

## Arguments

- moo_counts:

  counts dataframe or `multiOmicDataSet` containing `count_type` &
  `sub_count_type` in the counts slot

- count_type:

  the type of counts to use when `moo_counts` is a `multiOmicDataSet`;
  ignored for data frame input.

- sub_count_type:

  used when `count_type` refers to a list of count matrices; ignored for
  data frame input.

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
  you wish them to appear on the PCA plot. If `NULL`, no labels are
  added to PCA points. This can be the same Sample Names Column.
  However, you may desire different labels to display on your figure
  (e.g. shorter labels are sometimes preferred on plots). In that case,
  select the column with your preferred Labels here. The selected column
  should contain unique names for each sample.

- samples_to_rename:

  If you do not have a Plot Labels Column in your sample metadata table,
  you can use this parameter to rename samples manually for display on
  the PCA plot. Use "Add item" to add each additional sample for
  renaming. Use the following format to describe which old name (in your
  sample metadata table) you want to rename to which new name: old_name:
  new_name

- color_values:

  vector of colors as hex values or names recognized by R. Unnamed
  colors are assigned by factor level order when the grouping column is
  a factor; otherwise, they follow the order in which groups first
  appear in the metadata column. Defaults to `NULL`; when `NULL`,
  `mosuite_palette` is used for `data.frame` dispatch and stored colors
  are used for `multiOmicDataSet` dispatch.

- principal_components:

  vector with numbered principal components to plot

- legend_position:

  passed to in `legend.position`
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)

- point_size:

  size for
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)

- legend_font_size:

  font size for the PCA legend text. If `NULL`, the size is scaled
  automatically based on the number and length of legend labels.

- label_font_size:

  font size for text labels on the PCA plot.

- label_offset_x\_:

  horizontal offset for text labels on the PCA plot.

- label_offset_y\_:

  vertical offset for text labels on the PCA plot.

- log_transform:

  If `TRUE`, apply
  `log(x + log_transform_pseudocount, base = log_transform_base)` to
  sample count columns before PCA. Use this for count-like data such as
  raw, clean, filt, or CPM-like counts; leave it `FALSE` for already
  normalized/log-scale or batch-corrected values to avoid double
  transformation.

- log_transform_pseudocount:

  Pseudocount added before log-transforming counts when `log_transform`
  is `TRUE`.

- log_transform_base:

  Logarithm base to use when `log_transform` is `TRUE`. Use a numeric
  value, or `"e"`, `"ln"`, or `"natural"` for natural log. Default is
  `"ln"` to match the original PCA transform.

- interactive_plots:

  set to TRUE to make the PCA plot interactive with `plotly`.

- plots_subdir:

  subdirectory in `figures/` where PCA plots are saved.

- plot_filename:

  output filename for the PCA plot image.

- print_plots:

  Whether to print plots during analysis (Defaults to `FALSE`,
  overwritable using option 'moo_print_plots' or environment variable
  'MOO_PRINT_PLOTS')

- save_plots:

  Whether to save plots to files during analysis (Defaults to `TRUE`,
  overwritable using option 'moo_save_plots' or environment variable
  'MOO_SAVE_PLOTS')

- ...:

  arguments forwarded to method

## Value

ggplot object

## See also

[`plot_pca()`](https://ccbr.github.io/MOSuite/reference/plot_pca.md)
generic

Other PCA functions:
[`calc_pca()`](https://ccbr.github.io/MOSuite/reference/calc_pca.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/reference/plot_pca.md),
[`plot_pca_3d()`](https://ccbr.github.io/MOSuite/reference/plot_pca_3d.md)
