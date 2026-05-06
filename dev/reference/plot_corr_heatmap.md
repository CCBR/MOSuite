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
      color_values = c(
        "#5954d6", "#e1562c", "#b80058", "#00c6f8", "#d163e6", "#00a76c",
        "#ff9287", "#008cf9", "#006e00", "#796880", "#FFA500", "#878500"
      ))

## See also

- [`plot_corr_heatmap.multiOmicDataSet()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap-multiOmicDataSet.md)

- [`plot_corr_heatmap.data.frame()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap-data.frame.md)

Other plotters:
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap.md),
[`plot_histogram()`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.md),
[`print_or_save_plot()`](https://ccbr.github.io/MOSuite/dev/reference/print_or_save_plot.md)

Other heatmaps:
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap.md)

Other moo methods:
[`batch_correct_counts()`](https://ccbr.github.io/MOSuite/dev/reference/batch_correct_counts.md),
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/dev/reference/clean_raw_counts.md),
[`diff_counts()`](https://ccbr.github.io/MOSuite/dev/reference/diff_counts.md),
[`filter_counts()`](https://ccbr.github.io/MOSuite/dev/reference/filter_counts.md),
[`filter_diff()`](https://ccbr.github.io/MOSuite/dev/reference/filter_diff.md),
[`normalize_counts()`](https://ccbr.github.io/MOSuite/dev/reference/normalize_counts.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap.md),
[`plot_histogram()`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.md),
[`run_deseq2()`](https://ccbr.github.io/MOSuite/dev/reference/run_deseq2.md),
[`set_color_pal()`](https://ccbr.github.io/MOSuite/dev/reference/set_color_pal.md)

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
