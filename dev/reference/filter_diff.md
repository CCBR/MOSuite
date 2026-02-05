# Filter features from differential analysis based on statistical significance

Outputs dataset of significant genes from DEG table; filters genes based
on statistical significance (p-value or adjusted p-value) and change
(fold change, log2 fold change, or t-statistic); in addition allows for
selection of DEG estimates and for sub-setting of contrasts and groups
included in the output gene list.

## Usage

``` r
filter_diff(
  moo,
  feature_id_colname = NULL,
  significance_column = "adjpval",
  significance_cutoff = 0.05,
  change_column = "logFC",
  change_cutoff = 1,
  filtering_mode = "any",
  include_estimates = c("FC", "logFC", "tstat", "pval", "adjpval"),
  round_estimates = TRUE,
  rounding_decimal_for_percent_cells = 0,
  contrast_filter = "none",
  contrasts = c(),
  groups = c(),
  groups_filter = "none",
  label_font_size = 6,
  label_distance = 1,
  y_axis_expansion = 0.08,
  fill_colors = c("steelblue1", "whitesmoke"),
  pie_chart_in_3d = TRUE,
  bar_width = 0.4,
  draw_bar_border = TRUE,
  plot_type = "bar",
  plot_titles_fontsize = 12,
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = file.path("diff", "filt")
)
```

## Arguments

- moo:

  multiOmicDataSet object (see
  [`create_multiOmicDataSet_from_dataframes()`](https://ccbr.github.io/MOSuite/dev/reference/create_multiOmicDataSet_from_dataframes.md))

- feature_id_colname:

  The column from the counts data containing the Feature IDs (Usually
  Gene or Protein ID). This is usually the first column of your input
  Counts Matrix. Only columns of Text type from your input Counts Matrix
  will be available to select for this parameter. (Default: `NULL` -
  first column in the counts matrix will be used.)

- significance_column:

  Column name for significance, e.g. `"pval"` or `"pvaladj"` (default)

- significance_cutoff:

  Features will only be kept if their `significance_column` is less then
  this cutoff threshold

- change_column:

  Column name for change, e.g. `"logFC"` (default)

- change_cutoff:

  Features will only be kept if the absolute value of their
  `change_column` is greater than or equal to this cutoff threshold

- filtering_mode:

  Accepted values: `"any"` or `"all"` to include features that meet the
  criteria in *any* contrast or in *all* contrasts

- include_estimates:

  Column names of estimates to include. Default:
  `c("FC", "logFC", "tstat", "pval", "adjpval")`

- round_estimates:

  Whether to round estimates. Default: `TRUE`

- rounding_decimal_for_percent_cells:

  Decimal place to use when rounding Percent cells

- contrast_filter:

  Whether to filter `contrasts` in or our of analysis. If `"keep"`, only
  the contrast names listed in `contrasts` will be included. If
  `"remove`, the contrast names listed by `contrasts` will be removed.
  If `"none"`, all contrasts in the dataset are used. Options: `"keep"`,
  `"remove"`, or `"none"`

- contrasts:

  Contrast names to filter by `contrast_filter`. If `contrast_filter` is
  `"none"`, this parameter has no effect.

- groups:

  Group names to filter by `groups_filter`. If `groups_filter` is
  `"none"`, this parameter has no effect. Options: `"keep"`, `"remove"`,
  or `"none"`

- groups_filter:

  Whether to filter `groups` in or out of analysis. If `"keep"`, only
  the group names listed in `groups` will be included. If `"remove"`,
  the group names listed by `groups` will be removed. If `"none"`, all
  groups in the dataset are used.

- label_font_size:

  Font size for labels in the plot (default: 6)

- label_distance:

  Distance of labels from the bars (default: 1)

- y_axis_expansion:

  Expansion of the y-axis (default: 0.08)

- fill_colors:

  Fill colors for the bars (default: c("steelblue1", "whitesmoke"))

- pie_chart_in_3d:

  Whether to draw pie charts in 3D (default: TRUE)

- bar_width:

  Width of the bars (default: 0.4)

- draw_bar_border:

  Whether to draw borders around bars (default: TRUE)

- plot_type:

  "bar" or "pie"

- plot_titles_fontsize:

  Font size for plot titles (default: 12)

- print_plots:

  Whether to print plots during analysis (Defaults to `FALSE`,
  overwritable using option 'moo_print_plots' or environment variable
  'MOO_PRINT_PLOTS')

- save_plots:

  Whether to save plots to files during analysis (Defaults to `TRUE`,
  overwritable using option 'moo_save_plots' or environment variable
  'MOO_SAVE_PLOTS')

- plots_subdir:

  subdirectory in where plots will be saved if `save_plots` is `TRUE`

## See also

Other moo methods:
[`batch_correct_counts()`](https://ccbr.github.io/MOSuite/dev/reference/batch_correct_counts.md),
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/dev/reference/clean_raw_counts.md),
[`diff_counts()`](https://ccbr.github.io/MOSuite/dev/reference/diff_counts.md),
[`filter_counts()`](https://ccbr.github.io/MOSuite/dev/reference/filter_counts.md),
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
    "filt" = as.data.frame(nidap_filtered_counts)
  )
) |>
  diff_counts(
    count_type = "filt",
    sub_count_type = NULL,
    sample_id_colname = "Sample",
    feature_id_colname = "Gene",
    covariates_colnames = c("Group", "Batch"),
    contrast_colname = c("Group"),
    contrasts = c("B-A", "C-A", "B-C"),
    voom_normalization_method = "quantile",
  ) |>
  filter_diff()
#> * differential counts
#> Setting first column of `counts` as gene annotation.
#> Total number of genes included: 7943
#> Saving 6.67 x 6.67 in image
#> `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'
#> Joining with `by = join_by(Gene)`
#> Joining with `by = join_by(Gene)`
#> * filtering differential features
#> Total number of genes selected with adjpval < 0.05 and | logFC | ≥ 1 is sum(selgenes)
#> Saving 6.67 x 6.67 in image
head(moo@analyses$diff_filt)
#>            Gene B-A_FC B-A_logFC B-A_tstat B-A_pval B-A_adjpval C-A_FC
#> 1 1110034G24Rik  21.70     4.440      3.20  0.00782       0.210  36.60
#> 2 3110082I17Rik  -1.73    -0.789     -1.35  0.20300       0.710 -21.90
#> 3 4632428N05Rik   2.43     1.280      2.76  0.01770       0.303   4.66
#> 4 4833439L19Rik  -1.38    -0.460     -1.18  0.26000       0.758  -3.59
#> 5 4930523C07Rik  -2.30    -1.200     -1.62  0.13300       0.617   4.50
#> 6 5430427O19Rik  -2.22    -1.150     -2.46  0.03070       0.377  -4.49
#>   C-A_logFC C-A_tstat C-A_pval C-A_adjpval B-C_FC B-C_logFC B-C_tstat B-C_pval
#> 1      5.20      4.15 0.001410     0.02700  -1.69    -0.758    -0.838 0.419000
#> 2     -4.46     -3.80 0.002650     0.03830  12.70     3.670     2.930 0.012900
#> 3      2.22      5.25 0.000222     0.00929  -1.92    -0.941    -3.150 0.008590
#> 4     -1.84     -3.76 0.002810     0.03950   2.61     1.380     2.630 0.022200
#> 5      2.17      4.50 0.000767     0.01910 -10.30    -3.370    -5.040 0.000311
#> 6     -2.17     -3.68 0.003270     0.04320   2.02     1.010     1.530 0.153000
#>   B-C_adjpval
#> 1      0.7070
#> 2      0.1440
#> 3      0.1240
#> 4      0.1860
#> 5      0.0224
#> 6      0.4420
```
