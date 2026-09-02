# Extract count data

Re-exported from
[MOObject::extract_counts](https://ccbr.github.io/MOObject/reference/extract_counts.html)
– view the MOObject docs for details.

## Usage

``` r
extract_counts(moo, count_type, sub_count_type = NULL)
```

## Arguments

- moo:

  multiOmicDataSet containing `count_type` & `sub_count_type` in the
  counts slot

- count_type:

  the type of counts to use – must be a name in the counts slot
  (`moo@counts[[count_type]]`)

- sub_count_type:

  if `count_type` is a list, specify the sub count type within the list
  (`moo@counts[[count_type]][[sub_count_type]]`). (Default: `NULL`)

## Value

A data frame of counts.

## See also

[`MOObject::extract_counts()`](https://ccbr.github.io/MOObject/reference/extract_counts.html)

Other moo methods:
[`batch_correct_counts()`](https://ccbr.github.io/MOSuite/dev/reference/batch_correct_counts.md),
[`calc_cpm()`](https://ccbr.github.io/MOSuite/dev/reference/calc_cpm.md),
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/dev/reference/clean_raw_counts.md),
[`diff_counts()`](https://ccbr.github.io/MOSuite/dev/reference/diff_counts.md),
[`filter_counts()`](https://ccbr.github.io/MOSuite/dev/reference/filter_counts.md),
[`filter_diff()`](https://ccbr.github.io/MOSuite/dev/reference/filter_diff.md),
[`normalize_counts()`](https://ccbr.github.io/MOSuite/dev/reference/normalize_counts.md),
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap.md),
[`plot_histogram()`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md),
[`plot_pca_2d()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_2d.md),
[`plot_pca_3d()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_3d.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.md),
[`plot_venn_diagram()`](https://ccbr.github.io/MOSuite/dev/reference/plot_venn_diagram.md),
[`plot_volcano_enhanced()`](https://ccbr.github.io/MOSuite/dev/reference/plot_volcano_enhanced.md),
[`plot_volcano_summary()`](https://ccbr.github.io/MOSuite/dev/reference/plot_volcano_summary.md),
[`run_deseq2()`](https://ccbr.github.io/MOSuite/dev/reference/run_deseq2.md),
[`set_color_pal()`](https://ccbr.github.io/MOSuite/dev/reference/set_color_pal.md),
[`set_default_colors()`](https://ccbr.github.io/MOSuite/dev/reference/set_default_colors.md)
