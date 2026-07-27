# Perform batch correction

Perform batch correction using sva::ComBat()

## Usage

``` r
batch_correct_counts(
  moo,
  count_type = "norm",
  sub_count_type = "voom",
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  samples_to_include = NULL,
  covariates_colnames = "Group",
  batch_colname = "Batch",
  label_colname = "Label",
  samples_to_rename = c(""),
  add_label_to_pca = TRUE,
  principal_component_on_x_axis = 1,
  principal_component_on_y_axis = 2,
  legend_position_for_pca = "top",
  label_offset_x_ = 2,
  label_offset_y_ = 2,
  label_font_size = 3,
  point_size_for_pca = 5,
  color_histogram_by_group = TRUE,
  set_min_max_for_x_axis_for_histogram = FALSE,
  minimum_for_x_axis_for_histogram = -1,
  maximum_for_x_axis_for_histogram = 1,
  legend_font_size_for_histogram = NULL,
  legend_position_for_histogram = "top",
  number_of_histogram_legend_columns = 6,
  plot_corr_matrix_heatmap = TRUE,
  colors_for_plots = NULL,
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  interactive_plots = FALSE,
  plots_subdir = "batch"
)
```

## Arguments

- moo:

  multiOmicDataSet object (see
  [`create_multiOmicDataSet_from_dataframes()`](https://ccbr.github.io/MOSuite/reference/create_multiOmicDataSet_from_dataframes.md))

- count_type:

  the type of counts to use – must be a name in the counts slot
  (`moo@counts`)

- sub_count_type:

  if `count_type` is a list, specify the sub count type within the list.
  (Default: `"voom"`)

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

- samples_to_include:

  Which samples would you like to include? Usually, you will choose all
  sample columns, or you could choose to remove certain samples. Samples
  excluded here will be removed in this step and from further analysis
  downstream of this step. (Default: `NULL` - all sample IDs in
  `moo@sample_meta` will be used.)

- covariates_colnames:

  The column name(s) from the sample metadata containing variable(s) of
  interest, such as phenotype. Most commonly this will be the same
  column selected for your Groups Column. Some experimental designs may
  require that you add additional covariate columns here. Do not include
  the `batch_colname` here.

- batch_colname:

  The column from the sample metadata containing the batch information.
  Samples extracted, prepared, or sequenced at separate times or using
  separate materials/staff/equipment may belong to different batches.
  Not all data sets have batches, in which case you do not need batch
  correction. If your data set has no batches, you can provide a batch
  column with the same value in every row to skip batch correction
  (alternatively, simply do not run this function).

- label_colname:

  The column from the sample metadata containing the sample labels as
  you wish them to appear in heatmap and PCA figures. This can be the
  same Sample Names Column. However, you may desire different labels to
  display on your figures (e.g. shorter labels are sometimes preferred
  on plots). In that case, select the column with your preferred Labels
  here. The selected column should contain unique names for each sample.
  Use `add_label_to_pca` to control whether these labels are displayed
  on the PCA plot.

- samples_to_rename:

  If you do not have a Plot Labels Column in your sample metadata table,
  you can use this parameter to rename samples manually for display on
  the PCA plot. Use "Add item" to add each additional sample for
  renaming. Use the following format to describe which old name (in your
  sample metadata table) you want to rename to which new name: old_name:
  new_name

- add_label_to_pca:

  If `TRUE`, display labels from `label_colname` on PCA points. If
  `FALSE`, the PCA plot uses unlabeled points while heatmap labels still
  use `label_colname`.

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
  Color Histogram" parameter (below). Default is TRUE.

- set_min_max_for_x_axis_for_histogram:

  whether to set min/max value for histogram x-axis

- minimum_for_x_axis_for_histogram:

  x-axis minimum for histogram plot

- maximum_for_x_axis_for_histogram:

  x-axis maximum for histogram plot

- legend_font_size_for_histogram:

  legend font size for the histogram plot. If `NULL`, the size is scaled
  automatically.

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

  Optional colors for PCA/histogram/heatmap plots. If `NULL`, colors are
  taken from `moo@analyses$colors[[group_colname]]`. Colors must either
  be names in
  [`grDevices::colors()`](https://rdrr.io/r/grDevices/colors.html) or
  valid hex codes. Unnamed colors are assigned by factor level order
  when the grouping column is a factor; otherwise, they follow the order
  in which groups first appear in the metadata column. If more groups
  are present than colors provided, supplied colors are used first and
  additional colors are generated from the selected palette for the
  remaining groups; random colors are used only if that palette returns
  fewer colors than the number of groups.

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

`multiOmicDataSet` with batch-corrected counts

## See also

Other moo methods:
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/reference/clean_raw_counts.md),
[`diff_counts()`](https://ccbr.github.io/MOSuite/reference/diff_counts.md),
[`filter_counts()`](https://ccbr.github.io/MOSuite/reference/filter_counts.md),
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
moo <- multiOmicDataSet(
  sample_metadata = as.data.frame(nidap_sample_metadata),
  anno_dat = data.frame(),
  counts_lst = list(
    "raw" = as.data.frame(nidap_raw_counts),
    "clean" = as.data.frame(nidap_clean_raw_counts),
    "filt" = as.data.frame(nidap_filtered_counts),
    "norm" = list(
      "voom" = as.data.frame(nidap_norm_counts)
    )
  )
) |>
  batch_correct_counts(
    count_type = "norm",
    sub_count_type = "voom",
    covariates_colnames = "Group",
    batch_colname = "Batch",
    label_colname = "Label"
  )
#> * batch-correcting norm-voom counts
#> Found2batches
#> Adjusting for2covariate(s) or covariate level(s)
#> Standardizing Data across genes
#> Fitting L/S model and finding priors
#> Finding parametric adjustments
#> Adjusting the Data
#> Saving 6.67 x 6.67 in image
#> The total number of features in output: 7943
#> Number of samples after batch correction: 10

head(moo@counts[["batch"]])
#>            Gene       A1       A2       A3       B1       B2       B3       C1
#> 1 0610007P14Rik 6.437738 6.251229 6.048600 6.284429 6.188062 6.180803 6.333751
#> 2 0610009B22Rik 4.904608 5.100317 4.960486 4.037742 4.843373 5.098318 4.013808
#> 3 0610010F05Rik 4.921026 5.701279 6.485933 6.140332 5.847360 5.560233 3.737422
#> 4 0610011F06Rik 5.309874 5.288411 5.069086 5.261067 5.269024 5.551350 5.548404
#> 5 0610012G03Rik 5.426686 5.406358 5.415468 4.625768 5.333482 5.529869 5.845995
#> 6 0610037L13Rik 5.413417 5.293344 5.144240 5.421276 3.945936 4.831507 4.443280
#>         C2       C3
#> 1 6.253867 6.530433
#> 2 4.391701 5.050022
#> 3 2.756696 2.865261
#> 4 5.919472 5.455400
#> 5 6.086350 4.769502
#> 6 4.651311 5.063511
```
