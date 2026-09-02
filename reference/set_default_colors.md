# Set default colors on a multiOmicDataSet

Populates `moo@analyses$colors` using
[`get_colors_lst()`](https://ccbr.github.io/MOSuite/reference/get_colors_lst.md)
with the default MOSuite palette.

## Usage

``` r
set_default_colors(moo, palette = mosuite_palette)
```

## Arguments

- moo:

  `multiOmicDataSet` object.

- palette:

  Character vector of colors to assign. Defaults to `mosuite_palette`.

## Value

`moo` with constructor-style default colors set at
`moo@analyses$colors`.

## See also

Other moo methods:
[`batch_correct_counts()`](https://ccbr.github.io/MOSuite/reference/batch_correct_counts.md),
[`calc_cpm()`](https://ccbr.github.io/MOSuite/reference/calc_cpm.md),
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/reference/clean_raw_counts.md),
[`diff_counts()`](https://ccbr.github.io/MOSuite/reference/diff_counts.md),
[`extract_counts()`](https://ccbr.github.io/MOSuite/reference/extract_counts.md),
[`filter_counts()`](https://ccbr.github.io/MOSuite/reference/filter_counts.md),
[`filter_diff()`](https://ccbr.github.io/MOSuite/reference/filter_diff.md),
[`normalize_counts()`](https://ccbr.github.io/MOSuite/reference/normalize_counts.md),
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_corr_heatmap.md),
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
[`set_color_pal()`](https://ccbr.github.io/MOSuite/reference/set_color_pal.md)

## Examples

``` r
moo <- multiOmicDataSet(
  sample_metadata = as.data.frame(nidap_sample_metadata),
  anno_dat = data.frame(),
  counts_lst = list("raw" = as.data.frame(nidap_raw_counts))
)
moo <- set_default_colors(moo)
names(moo@analyses$colors)
#> [1] "Sample"    "Group"     "Replicate" "Batch"     "Label"    
```
