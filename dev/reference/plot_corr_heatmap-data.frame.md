# Plot correlation heatmap for counts dataframe

Plot correlation heatmap for counts dataframe

## Arguments

- moo_counts:

  a `data.frame` of counts

- sample_metadata:

  sample metadata as a data frame or tibble (**Required**)

- sample_id_colname:

  The column from the sample metadata containing the sample names. The
  names in this column must exactly match the names used as the sample
  column names of your input Counts Matrix. (Default: `NULL` - first
  column in the sample metadata will be used.)

- feature_id_colname:

  The column from the counts data containing the Feature IDs (Usually
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

## See also

[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap.md)
generic

Other plotters for counts dataframes:
[`plot_histogram,data.frame-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram.data.frame.md),
[`plot_pca,data.frame-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.data.frame.md),
[`plot_read_depth,data.frame-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.data.frame.md)
