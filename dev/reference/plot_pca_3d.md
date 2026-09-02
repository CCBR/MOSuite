# Perform and plot a 3D Principal Components Analysis

Perform and plot a 3D Principal Components Analysis

3D PCA for counts dataframe

## Usage

``` r
plot_pca_3d(
  moo_counts,
  count_type = NULL,
  sub_count_type = NULL,
  sample_metadata = NULL,
  feature_id_colname = NULL,
  sample_id_colname = NULL,
  samples_to_rename = NULL,
  group_colname = "Group",
  label_colname = "Label",
  principal_components = c(1, 2, 3),
  point_size = 8,
  label_font_size = 24,
  log_transform = FALSE,
  log_transform_pseudocount = 0.5,
  log_transform_base = "ln",
  color_values = NULL,
  plot_title = "PCA 3D",
  plot_filename = "pca_3D.html",
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = "pca",
  ...
)

## S7 method for class <MOObject::multiOmicDataSet>
plot_pca_3d(
  moo_counts,
  count_type = NULL,
  sub_count_type = NULL,
  sample_metadata = NULL,
  feature_id_colname = NULL,
  sample_id_colname = NULL,
  samples_to_rename = NULL,
  group_colname = "Group",
  label_colname = "Label",
  principal_components = c(1, 2, 3),
  point_size = 8,
  label_font_size = 24,
  log_transform = FALSE,
  log_transform_pseudocount = 0.5,
  log_transform_base = "ln",
  color_values = NULL,
  plot_title = "PCA 3D",
  plot_filename = "pca_3D.html",
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = "pca",
  ...
)

## S7 method for class <data.frame>
plot_pca_3d(
  moo_counts,
  count_type = NULL,
  sub_count_type = NULL,
  sample_metadata = NULL,
  feature_id_colname = NULL,
  sample_id_colname = NULL,
  samples_to_rename = NULL,
  group_colname = "Group",
  label_colname = "Label",
  principal_components = c(1, 2, 3),
  point_size = 8,
  label_font_size = 24,
  log_transform = FALSE,
  log_transform_pseudocount = 0.5,
  log_transform_base = "ln",
  color_values = NULL,
  plot_title = "PCA 3D",
  plot_filename = "pca_3D.html",
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = "pca",
  ...
)
```

## Arguments

- moo_counts:

  counts dataframe

- count_type:

  the type of counts to use. Ignored when `moo_counts` is already a
  dataframe.

- sub_count_type:

  used if `count_type` is a list in the counts slot: specify the sub
  count type within the list.

- sample_metadata:

  sample metadata as a data frame or tibble.

- feature_id_colname:

  The column from the counts data containing feature IDs. If `NULL`,
  first column is used.

- sample_id_colname:

  The column from sample metadata containing sample names. If `NULL`,
  first column is used.

- samples_to_rename:

  optional named mapping in `old_name: new_name` format for display
  labels.

- group_colname:

  The column from sample metadata containing sample group information.

- label_colname:

  The column from sample metadata containing sample labels.

- principal_components:

  vector with numbered principal components to plot

- point_size:

  size for
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)

- label_font_size:

  font size used for labels in the interactive figure.

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

- color_values:

  vector of colors as hex values or names recognized by R. Unnamed
  colors are assigned by factor level order when the grouping column is
  a factor; otherwise, they follow the order in which groups first
  appear in the metadata column. Defaults to `NULL`; when `NULL`,
  `mosuite_palette` is used for `data.frame` dispatch and stored colors
  are used for `multiOmicDataSet` dispatch.

- plot_title:

  title for the plot

- plot_filename:

  output filename when saving plots.

- print_plots:

  whether to print plot to the active graphics device.

- save_plots:

  whether to save plot to disk.

- plots_subdir:

  output subdirectory for saved plots.

- ...:

  additional arguments passed to methods

## Value

[`plotly::plot_ly`](https://rdrr.io/pkg/plotly/man/plot_ly.html) figure

## See also

Other moo methods:
[`batch_correct_counts()`](https://ccbr.github.io/MOSuite/dev/reference/batch_correct_counts.md),
[`calc_cpm()`](https://ccbr.github.io/MOSuite/dev/reference/calc_cpm.md),
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/dev/reference/clean_raw_counts.md),
[`diff_counts()`](https://ccbr.github.io/MOSuite/dev/reference/diff_counts.md),
[`extract_counts()`](https://ccbr.github.io/MOSuite/dev/reference/extract_counts.md),
[`filter_counts()`](https://ccbr.github.io/MOSuite/dev/reference/filter_counts.md),
[`filter_diff()`](https://ccbr.github.io/MOSuite/dev/reference/filter_diff.md),
[`normalize_counts()`](https://ccbr.github.io/MOSuite/dev/reference/normalize_counts.md),
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap.md),
[`plot_histogram()`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md),
[`plot_pca_2d()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_2d.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.md),
[`plot_venn_diagram()`](https://ccbr.github.io/MOSuite/dev/reference/plot_venn_diagram.md),
[`plot_volcano_enhanced()`](https://ccbr.github.io/MOSuite/dev/reference/plot_volcano_enhanced.md),
[`plot_volcano_summary()`](https://ccbr.github.io/MOSuite/dev/reference/plot_volcano_summary.md),
[`run_deseq2()`](https://ccbr.github.io/MOSuite/dev/reference/run_deseq2.md),
[`set_color_pal()`](https://ccbr.github.io/MOSuite/dev/reference/set_color_pal.md),
[`set_default_colors()`](https://ccbr.github.io/MOSuite/dev/reference/set_default_colors.md)

Other PCA functions:
[`calc_pca()`](https://ccbr.github.io/MOSuite/dev/reference/calc_pca.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md),
[`plot_pca_2d()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_2d.md)
