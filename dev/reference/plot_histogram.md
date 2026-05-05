# Plot histogram

Plot histogram

Plot histogram for multiOmicDataSet

Plot histogram for counts dataframe

## Usage

``` r
plot_histogram(moo_counts, ...)
```

## Arguments

- moo_counts:

  counts dataframe or `multiOmicDataSet` containing `count_type` &
  `sub_count_type` in the counts slot

- ...:

  additional arguments (ignored; accepted for compatibility with the moo
  dispatch)

- count_type:

  Required if `moo_counts` is a `multiOmicDataSet`: the type of counts
  to use – must be a name in the counts slot (`moo@counts`).

- sub_count_type:

  Used if `moo_counts` is a `multiOmicDataSet` AND if `count_type` is a
  list, specify the sub count type within the list

- sample_metadata:

  sample metadata as a data frame or tibble (**required**)

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

- color_values:

  vector of colors as hex values or names recognized by R

- color_by_group:

  Set to FALSE to label histogram by Sample Names, or set to TRUE to
  label histogram by the column you select in the "Group Column Used to
  Color Histogram" parameter (below). Default is FALSE.

- set_min_max_for_x_axis:

  whether to override the default for
  [`ggplot2::xlim()`](https://ggplot2.tidyverse.org/reference/lims.html)
  (default: `FALSE`)

- minimum_for_x_axis:

  value to override default `min` for
  [`ggplot2::xlim()`](https://ggplot2.tidyverse.org/reference/lims.html)

- maximum_for_x_axis:

  value to override default `max` for
  [`ggplot2::xlim()`](https://ggplot2.tidyverse.org/reference/lims.html)

- x_axis_label:

  text label for the x axis
  [`ggplot2::xlab()`](https://ggplot2.tidyverse.org/reference/labs.html)

- y_axis_label:

  text label for the y axis
  [`ggplot2::ylab()`](https://ggplot2.tidyverse.org/reference/labs.html)

- legend_position:

  passed to in `legend.position`
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)

- legend_font_size:

  passed to
  [`ggplot2::element_text()`](https://ggplot2.tidyverse.org/reference/element.html)
  via
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)

- number_of_legend_columns:

  passed to `ncol` in
  [`ggplot2::guide_legend()`](https://ggplot2.tidyverse.org/reference/guide_legend.html)

- interactive_plots:

  set to TRUE to make the plot interactive with `plotly`, allowing you
  to hover your mouse over a point or line to view sample information.
  The similarity heat map will not display if this toggle is set to
  TRUE. Default is FALSE.

## Value

ggplot object

## See also

`plot_histogram()` generic

Other plotters:
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.md),
[`print_or_save_plot()`](https://ccbr.github.io/MOSuite/dev/reference/print_or_save_plot.md)

Other moo methods:
[`batch_correct_counts()`](https://ccbr.github.io/MOSuite/dev/reference/batch_correct_counts.md),
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/dev/reference/clean_raw_counts.md),
[`diff_counts()`](https://ccbr.github.io/MOSuite/dev/reference/diff_counts.md),
[`filter_counts()`](https://ccbr.github.io/MOSuite/dev/reference/filter_counts.md),
[`filter_diff()`](https://ccbr.github.io/MOSuite/dev/reference/filter_diff.md),
[`normalize_counts()`](https://ccbr.github.io/MOSuite/dev/reference/normalize_counts.md),
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.md),
[`run_deseq2()`](https://ccbr.github.io/MOSuite/dev/reference/run_deseq2.md),
[`set_color_pal()`](https://ccbr.github.io/MOSuite/dev/reference/set_color_pal.md)

Other plotters for multiOmicDataSets:
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.md)

Other plotters for counts dataframes:
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.md)

## Examples

``` r
# plot histogram for a counts slot in a multiOmicDataset Object
moo <- multiOmicDataSet(
  sample_metadata = nidap_sample_metadata,
  anno_dat = data.frame(),
  counts_lst = list("raw" = nidap_raw_counts)
)
p <- plot_histogram(moo, count_type = "raw")

# customize the plot
plot_histogram(moo,
  count_type = "raw",
  group_colname = "Group", color_by_group = TRUE
)


# plot histogram for a counts dataframe directly
counts_dat <- moo@counts$raw
plot_histogram(
  counts_dat,
  sample_metadata = nidap_sample_metadata,
  sample_id_colname = "Sample",
  feature_id_colname = "GeneName",
  label_colname = "Label"
)


# plot histogram for a counts slot in a multiOmicDataset Object
moo <- multiOmicDataSet(
  sample_metadata = nidap_sample_metadata,
  anno_dat = data.frame(),
  counts_lst = list("raw" = nidap_raw_counts)
)
p <- plot_histogram(moo, count_type = "raw")

# customize the plot
plot_histogram(moo,
  count_type = "raw",
  group_colname = "Group", color_by_group = TRUE
)



# plot histogram for a counts dataframe directly
plot_histogram(
  nidap_clean_raw_counts,
  sample_metadata = nidap_sample_metadata,
  sample_id_colname = "Sample",
  feature_id_colname = "Gene",
  label_colname = "Label"
)


# customize the plot
plot_histogram(
  nidap_clean_raw_counts,
  sample_metadata = nidap_sample_metadata,
  sample_id_colname = "Sample",
  feature_id_colname = "Gene",
  group_colname = "Group",
  color_by_group = TRUE
)

```
