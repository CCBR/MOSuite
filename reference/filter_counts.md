# Filter low counts

This is often the first step in the QC portion of an analysis to filter
out features that have very low raw counts across most or all of your
samples.

## Usage

``` r
filter_counts(
  moo,
  count_type = "clean",
  feature_id_colname = NULL,
  sample_id_colname = NULL,
  group_colname = "Group",
  label_colname = NULL,
  samples_to_include = NULL,
  minimum_count_value_to_be_considered_nonzero = 8,
  minimum_number_of_samples_with_nonzero_counts_in_total = 7,
  minimum_number_of_samples_with_nonzero_counts_in_a_group = 3,
  use_cpm_counts_to_filter = TRUE,
  use_group_based_filtering = FALSE,
  principal_component_on_x_axis = 1,
  principal_component_on_y_axis = 2,
  legend_position_for_pca = "top",
  point_size_for_pca = 1,
  add_label_to_pca = TRUE,
  label_font_size = 3,
  label_offset_y_ = 2,
  label_offset_x_ = 2,
  samples_to_rename = c(""),
  color_histogram_by_group = FALSE,
  set_min_max_for_x_axis_for_histogram = FALSE,
  minimum_for_x_axis_for_histogram = -1,
  maximum_for_x_axis_for_histogram = 1,
  legend_position_for_histogram = "top",
  legend_font_size_for_histogram = 10,
  number_of_histogram_legend_columns = 6,
  colors_for_plots = NULL,
  plot_corr_matrix_heatmap = TRUE,
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  interactive_plots = FALSE,
  plots_subdir = "filt"
)
```

## Arguments

- moo:

  multiOmicDataSet object (see
  [`create_multiOmicDataSet_from_dataframes()`](https://ccbr.github.io/MOSuite/reference/create_multiOmicDataSet_from_dataframes.md))

- count_type:

  the type of counts to use – must be a name in the counts slot
  (`moo@counts`)

- feature_id_colname:

  The column from the counts data containing the Feature IDs (Usually
  Gene or Protein ID). This is usually the first column of your input
  Counts Matrix. Only columns of Text type from your input Counts Matrix
  will be available to select for this parameter. (Default: `NULL` -
  first column in the counts matrix will be used.)

- sample_id_colname:

  The column from the sample metadata containing the sample names. The
  names in this column must exactly match the names used as the sample
  column names of your input Counts Matrix. (Default: `NULL` - first
  column in the sample metadata will be used.)

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

- samples_to_include:

  Which samples would you like to include? Usually, you will choose all
  sample columns, or you could choose to remove certain samples. Samples
  excluded here will be removed in this step and from further analysis
  downstream of this step. (Default: `NULL` - all sample IDs in
  `moo@sample_meta` will be used.)

- minimum_count_value_to_be_considered_nonzero:

  Minimum count value to be considered non-zero for a sample

- minimum_number_of_samples_with_nonzero_counts_in_total:

  Minimum number of samples (total) with non-zero counts

- minimum_number_of_samples_with_nonzero_counts_in_a_group:

  Only keeps genes that have at least this number of samples with
  nonzero CPM counts in at least one group

- use_cpm_counts_to_filter:

  If no transformation has been been performed on counts matrix (eg Raw
  Counts) set to TRUE. If TRUE counts will be transformed to CPM and
  filtered based on given criteria. If gene counts matrix has been
  transformed (eg log2, CPM, FPKM or some form of Normalization) set to
  FALSE. If FALSE no further transformation will be applied and features
  will be filtered as is. For RNAseq data RAW counts should be
  transformed to CPM in order to properly filter.

- use_group_based_filtering:

  If TRUE, only keeps features (e.g. genes) that have at least a certain
  number of samples with nonzero CPM counts in at least one group

- principal_component_on_x_axis:

  The principal component to plot on the x-axis for the PCA plot.
  Choices include 1, 2, 3, ... (default: 1)

- principal_component_on_y_axis:

  The principal component to plot on the y-axis for the PCA plot.
  Choices include 1, 2, 3, ... (default: 2)

- legend_position_for_pca:

  legend position for the PCA plot

- point_size_for_pca:

  geom point size for the PCA plot

- add_label_to_pca:

  label points on the PCA plot

- label_font_size:

  label font size for the PCA plot

- label_offset_y\_:

  label offset y for the PCA plot

- label_offset_x\_:

  label offset x for the PCA plot

- samples_to_rename:

  If you do not have a Plot Labels Column in your sample metadata table,
  you can use this parameter to rename samples manually for display on
  the PCA plot. Use "Add item" to add each additional sample for
  renaming. Use the following format to describe which old name (in your
  sample metadata table) you want to rename to which new name: old_name:
  new_name

- color_histogram_by_group:

  Set to FALSE to label histogram by Sample Names, or set to TRUE to
  label histogram by the column you select in the "Group Column Used to
  Color Histogram" parameter (below). Default is FALSE.

- set_min_max_for_x_axis_for_histogram:

  whether to set min/max value for histogram x-axis

- minimum_for_x_axis_for_histogram:

  x-axis minimum for histogram plot

- maximum_for_x_axis_for_histogram:

  x-axis maximum for histogram plot

- legend_position_for_histogram:

  legend position for the histogram plot. consider setting to 'none' for
  a large number of samples.

- legend_font_size_for_histogram:

  legend font size for the histogram plot

- number_of_histogram_legend_columns:

  number of columns for the histogram legend

- colors_for_plots:

  Colors for the PCA and histogram will be picked, in order, from this
  list. Colors must either be names in
  [`grDevices::colors()`](https://rdrr.io/r/grDevices/colors.html) or
  valid hex codes.

- plot_corr_matrix_heatmap:

  Datasets with a large number of samples may be too large to create a
  correlation matrix heatmap. If this function takes longer than 5
  minutes to run, Set to `FALSE` and the correlation matrix will not be
  be created. Default is `TRUE`.

- print_plots:

  Whether to print plots during analysis (Defaults to `FALSE`,
  overwritable using option 'moo_print_plots' or environment variable
  'MOO_PRINT_PLOTS')

- save_plots:

  Whether to save plots to files during analysis (Defaults to `TRUE`,
  overwritable using option 'moo_save_plots' or environment variable
  'MOO_SAVE_PLOTS')

- interactive_plots:

  set to TRUE to make PCA and Histogram plots interactive with `plotly`,
  allowing you to hover your mouse over a point or line to view sample
  information. The similarity heat map will not display if this toggle
  is set to `TRUE`. Default is `FALSE`.

- plots_subdir:

  subdirectory in `figures/` where plots will be saved if `save_plots`
  is `TRUE`

## Value

`multiOmicDataSet` with filtered counts

## Details

This function takes a multiOmicDataSet containing clean raw counts and a
sample metadata table, and returns the multiOmicDataSet object with
filtered counts. It also produces an image consisting of three QC plots.

You can tune the threshold for tuning how low counts for a given gene
are before they are deemed "too low" and filtered out of downstream
analysis. By default, this parameter is set to 1, meaning any raw count
value less than 1 will count as "too low".

The QC plots are provided to help you assess: (1) PCA Plot: the within
and between group variance in expression after dimensionality reduction;
(2) Count Density Histogram: the dis/similarity of count distributions
between samples; and (3) Similarity Heatmap: the overall similarity of
samples to one another based on unsupervised clustering.

## See also

Other moo methods:
[`batch_correct_counts()`](https://ccbr.github.io/MOSuite/reference/batch_correct_counts.md),
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/reference/clean_raw_counts.md),
[`diff_counts()`](https://ccbr.github.io/MOSuite/reference/diff_counts.md),
[`filter_diff()`](https://ccbr.github.io/MOSuite/reference/filter_diff.md),
[`normalize_counts()`](https://ccbr.github.io/MOSuite/reference/normalize_counts.md),
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_corr_heatmap.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_expr_heatmap.md),
[`plot_histogram()`](https://ccbr.github.io/MOSuite/reference/plot_histogram.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/reference/plot_pca.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/reference/plot_read_depth.md),
[`run_deseq2()`](https://ccbr.github.io/MOSuite/reference/run_deseq2.md),
[`set_color_pal()`](https://ccbr.github.io/MOSuite/reference/set_color_pal.md)

## Examples

``` r
moo <- create_multiOmicDataSet_from_dataframes(
  as.data.frame(nidap_sample_metadata),
  as.data.frame(nidap_clean_raw_counts),
  sample_id_colname = "Sample",
  feature_id_colname = "Gene"
) %>%
  filter_counts(
    count_type = "raw"
  )
#> * filtering raw counts
#> Number of features after filtering: 7943
#> colors_for_plots NULL
#> colors_for_plots character
#> Saving 6.67 x 6.67 in image
#> Saving 6.67 x 6.67 in image
head(moo@counts$filt)
#>      Gene   A1   A2   A3   B1   B2   B3   C1   C2   C3
#> 1  Mrpl15 1245 1341 1476  965 1235 1784 1058 1732 1531
#> 2  Lypla1 1483 1410 1370 1146 1422 2624  991 1101 2352
#> 3   Tcea1 1381 2044 2051 2325 2386 1893 2391  916 2261
#> 4 Atp6v1h 1033 1959 1890 2075 2702 2150 2436 1321 1018
#> 5  Rb1cc1  666 1397 1576  681 2040 1988  774 1921 2660
#> 6  Pcmtd1  798  966  407  487  455  950 1710 1995 2502
```
