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
  label_colname = NULL,
  colors_for_plots = NULL,
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = "batch"
)
```

## Arguments

- moo:

  multiOmicDataSet object (see
  [`create_multiOmicDataSet_from_dataframes()`](https://ccbr.github.io/MOSuite/dev/reference/create_multiOmicDataSet_from_dataframes.md))

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
  you wish them to appear in the plots produced by this template. This
  can be the same Sample Names Column. However, you may desire different
  labels to display on your figure (e.g. shorter labels are sometimes
  preferred on plots). In that case, select the column with your
  preferred Labels here. The selected column should contain unique names
  for each sample. (Default: `NULL` – `sample_id_colname` will be used.)

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

- plots_subdir:

  subdirectory in `figures/` where plots will be saved if `save_plots`
  is `TRUE`

## Value

`multiOmicDataSet` with batch-corrected counts

## See also

Other moo methods:
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/dev/reference/clean_raw_counts.md),
[`diff_counts()`](https://ccbr.github.io/MOSuite/dev/reference/diff_counts.md),
[`filter_counts()`](https://ccbr.github.io/MOSuite/dev/reference/filter_counts.md),
[`filter_diff()`](https://ccbr.github.io/MOSuite/dev/reference/filter_diff.md),
[`normalize_counts()`](https://ccbr.github.io/MOSuite/dev/reference/normalize_counts.md),
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap.md),
[`plot_histogram()`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.md),
[`run_deseq2()`](https://ccbr.github.io/MOSuite/dev/reference/run_deseq2.md),
[`set_color_pal()`](https://ccbr.github.io/MOSuite/dev/reference/set_color_pal.md)

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
) %>%
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
#> Saving 6.67 x 6.67 in image
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
