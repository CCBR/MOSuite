# Plot correlation heatmap

Plot correlation heatmap

## Usage

``` r
plot_corr_heatmap(moo_counts, ...)
```

## Arguments

- moo_counts:

  counts dataframe or `multiOmicDataSet` containing `count_type` &
  `sub_count_type` in the counts slot

- ...:

  arguments forwarded to method

## Value

heatmap from
[`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)

## Details

### Method Usage

    # multiOmicDataSet
    plot_corr_heatmap(moo_counts,
      count_type,
      sub_count_type = NULL,
      ...)

    # dataframe
    plot_corr_heatmap(moo_counts,
      sample_metadata,
      sample_id_colname = NULL,
      feature_id_colname = NULL,
      group_colname = "Group",
      label_colname = "Label",
      color_values = mosuite_palette)

## See also

- [`plot_corr_heatmap.multiOmicDataSet()`](https://ccbr.github.io/MOSuite/reference/plot_corr_heatmap-multiOmicDataSet.md)

- [`plot_corr_heatmap.data.frame()`](https://ccbr.github.io/MOSuite/reference/plot_corr_heatmap-data.frame.md)

Other plotters:
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_expr_heatmap.md),
[`plot_histogram()`](https://ccbr.github.io/MOSuite/reference/plot_histogram.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/reference/plot_pca.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/reference/plot_read_depth.md),
[`print_or_save_plot()`](https://ccbr.github.io/MOSuite/reference/print_or_save_plot.md)

Other heatmaps:
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_expr_heatmap.md)

Other moo methods:
[`batch_correct_counts()`](https://ccbr.github.io/MOSuite/reference/batch_correct_counts.md),
[`calc_cpm()`](https://ccbr.github.io/MOSuite/reference/calc_cpm.md),
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/reference/clean_raw_counts.md),
[`diff_counts()`](https://ccbr.github.io/MOSuite/reference/diff_counts.md),
[`extract_counts()`](https://ccbr.github.io/MOSuite/reference/extract_counts.md),
[`filter_counts()`](https://ccbr.github.io/MOSuite/reference/filter_counts.md),
[`filter_diff()`](https://ccbr.github.io/MOSuite/reference/filter_diff.md),
[`normalize_counts()`](https://ccbr.github.io/MOSuite/reference/normalize_counts.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_expr_heatmap.md),
[`plot_histogram()`](https://ccbr.github.io/MOSuite/reference/plot_histogram.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/reference/plot_pca.md),
[`plot_pca_2d()`](https://ccbr.github.io/MOSuite/reference/plot_pca_2d.md),
[`plot_pca_3d()`](https://ccbr.github.io/MOSuite/reference/plot_pca_3d.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/reference/plot_read_depth.md),
[`plot_venn_diagram()`](https://ccbr.github.io/MOSuite/reference/plot_venn_diagram.md),
[`plot_volcano_enhanced()`](https://ccbr.github.io/MOSuite/reference/plot_volcano_enhanced.md),
[`plot_volcano_summary()`](https://ccbr.github.io/MOSuite/reference/plot_volcano_summary.md),
[`run_deseq2()`](https://ccbr.github.io/MOSuite/reference/run_deseq2.md),
[`set_color_pal()`](https://ccbr.github.io/MOSuite/reference/set_color_pal.md),
[`set_default_colors()`](https://ccbr.github.io/MOSuite/reference/set_default_colors.md)

## Examples

``` r
# plot correlation heatmap for a counts slot in a multiOmicDataset Object
moo <- multiOmicDataSet(
  sample_metadata = as.data.frame(nidap_sample_metadata),
  anno_dat = data.frame(),
  counts_lst = list("raw" = as.data.frame(nidap_raw_counts))
)
p <- plot_corr_heatmap(moo, count_type = "raw")

# plot correlation heatmap for a counts dataframe
plot_corr_heatmap(
  moo@counts$raw,
  sample_metadata = moo@sample_meta,
  sample_id_colname = "Sample",
  feature_id_colname = "Gene",
  group_colname = "Group",
  label_colname = "Label"
)
```
