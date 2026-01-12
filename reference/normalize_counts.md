# Normalize counts

Normalize counts

## Usage

``` r
normalize_counts(
  moo,
  count_type = "filt",
  norm_type = "voom",
  feature_id_colname = NULL,
  samples_to_include = NULL,
  sample_id_colname = NULL,
  group_colname = "Group",
  label_colname = NULL,
  input_in_log_counts = FALSE,
  voom_normalization_method = "quantile",
  samples_to_rename = c(""),
  add_label_to_pca = TRUE,
  principal_component_on_x_axis = 1,
  principal_component_on_y_axis = 2,
  legend_position_for_pca = "top",
  label_offset_x_ = 2,
  label_offset_y_ = 2,
  label_font_size = 3,
  point_size_for_pca = 8,
  color_histogram_by_group = TRUE,
  set_min_max_for_x_axis_for_histogram = FALSE,
  minimum_for_x_axis_for_histogram = -1,
  maximum_for_x_axis_for_histogram = 1,
  legend_font_size_for_histogram = 10,
  legend_position_for_histogram = "top",
  number_of_histogram_legend_columns = 6,
  plot_corr_matrix_heatmap = TRUE,
  colors_for_plots = NULL,
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  interactive_plots = FALSE,
  plots_subdir = "norm"
)
```

## Arguments

- moo:

  multiOmicDataSet object (see
  [`create_multiOmicDataSet_from_dataframes()`](https://ccbr.github.io/MOSuite/reference/create_multiOmicDataSet_from_dataframes.md))

- count_type:

  the type of counts to use – must be a name in the counts slot
  (`moo@counts`)

- norm_type:

  normalization type. Default: "voom" which uses
  [`limma::voom`](https://rdrr.io/pkg/limma/man/voom.html).

- feature_id_colname:

  The column from the counts data containing the Feature IDs (Usually
  Gene or Protein ID). This is usually the first column of your input
  Counts Matrix. Only columns of Text type from your input Counts Matrix
  will be available to select for this parameter. (Default: `NULL` -
  first column in the counts matrix will be used.)

- samples_to_include:

  Which samples would you like to include? Usually, you will choose all
  sample columns, or you could choose to remove certain samples. Samples
  excluded here will be removed in this step and from further analysis
  downstream of this step. (Default: `NULL` - all sample IDs in
  `moo@sample_meta` will be used.)

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

- input_in_log_counts:

  set this to `TRUE` if counts are already log2-transformed

- voom_normalization_method:

  Normalization method to be applied to the logCPM values when using
  [`limma::voom`](https://rdrr.io/pkg/limma/man/voom.html)

- samples_to_rename:

  If you do not have a Plot Labels Column in your sample metadata table,
  you can use this parameter to rename samples manually for display on
  the PCA plot. Use "Add item" to add each additional sample for
  renaming. Use the following format to describe which old name (in your
  sample metadata table) you want to rename to which new name: old_name:
  new_name

- add_label_to_pca:

  label points on the PCA plot

- principal_component_on_x_axis:

  The principal component to plot on the x-axis for the PCA plot.
  Choices include 1, 2, 3, ... (default: 1)

- principal_component_on_y_axis:

  The principal component to plot on the y-axis for the PCA plot.
  Choices include 1, 2, 3, ... (default: 2)

- legend_position_for_pca:

  legend position for the PCA plot

- label_offset_x\_:

  label offset x for the PCA plot

- label_offset_y\_:

  label offset y for the PCA plot

- label_font_size:

  label font size for the PCA plot

- point_size_for_pca:

  geom point size for the PCA plot

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

- legend_font_size_for_histogram:

  legend font size for the histogram plot

- legend_position_for_histogram:

  legend position for the histogram plot. consider setting to 'none' for
  a large number of samples.

- number_of_histogram_legend_columns:

  number of columns for the histogram legend

- plot_corr_matrix_heatmap:

  Datasets with a large number of samples may be too large to create a
  correlation matrix heatmap. If this function takes longer than 5
  minutes to run, Set to `FALSE` and the correlation matrix will not be
  be created. Default is `TRUE`.

- colors_for_plots:

  Colors for the PCA and histogram will be picked, in order, from this
  list. Colors must either be names in
  [`grDevices::colors()`](https://rdrr.io/r/grDevices/colors.html) or
  valid hex codes.

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

`multiOmicDataSet` with normalized counts

## See also

Other moo methods:
[`batch_correct_counts()`](https://ccbr.github.io/MOSuite/reference/batch_correct_counts.md),
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/reference/clean_raw_counts.md),
[`diff_counts()`](https://ccbr.github.io/MOSuite/reference/diff_counts.md),
[`filter_counts()`](https://ccbr.github.io/MOSuite/reference/filter_counts.md),
[`filter_diff()`](https://ccbr.github.io/MOSuite/reference/filter_diff.md),
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_corr_heatmap.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_expr_heatmap.md),
[`plot_histogram()`](https://ccbr.github.io/MOSuite/reference/plot_histogram.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/reference/plot_pca.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/reference/plot_read_depth.md),
[`run_deseq2()`](https://ccbr.github.io/MOSuite/reference/run_deseq2.md),
[`set_color_pal()`](https://ccbr.github.io/MOSuite/reference/set_color_pal.md)

## Examples

``` r
moo <- multiOmicDataSet(
  sample_metadata = as.data.frame(nidap_sample_metadata),
  anno_dat = data.frame(),
  counts_lst = list(
    "raw" = as.data.frame(nidap_raw_counts),
    "clean" = as.data.frame(nidap_clean_raw_counts),
    "filt" = as.data.frame(nidap_filtered_counts)
  )
) %>%
  normalize_counts(
    group_colname = "Group",
    label_colname = "Label"
  )
#> * normalizing filt counts
#> Total number of features included: 7943
#> Saving 6.67 x 6.67 in image
#> Saving 6.67 x 6.67 in image
#> Sample columns: A1, Sample columns: A2, Sample columns: A3, Sample columns: B1, Sample columns: B2, Sample columns: B3, Sample columns: C1, Sample columns: C2, Sample columns: C3
head(moo@counts[["norm"]][["voom"]])
#>            Gene       A1       A2       A3       B1       B2       B3       C1
#> 1 0610007P14Rik 6.532994 6.192871 5.954869 6.375896 6.275880 6.119449 6.419913
#> 2 0610009B22Rik 4.484983 5.448875 5.286875 3.445612 4.451347 5.473886 3.500359
#> 3 0610010F05Rik 4.883688 5.668494 6.537590 6.216408 5.893089 5.498884 3.845207
#> 4 0610011F06Rik 5.199684 5.374085 5.112952 5.155558 5.163359 5.650929 5.441965
#> 5 0610012G03Rik 5.368118 5.445918 5.456511 4.567138 5.274928 5.625039 5.787457
#> 6 0610037L13Rik 5.327987 5.388747 5.233520 5.450169 3.656585 4.929386 4.274944
#>         C2       C3
#> 1 6.172204 6.497050
#> 2 4.709254 5.471951
#> 3 2.685177 2.805426
#> 4 6.043492 5.490958
#> 5 6.214163 4.682896
#> 6 4.744405 5.173531
```
