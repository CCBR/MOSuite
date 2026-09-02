# Calculate counts-per-million (CPM) on raw counts in a multiOmicDataSet

Calculate counts-per-million (CPM) on raw counts in a multiOmicDataSet

## Usage

``` r
calc_cpm(moo, ...)
```

## Arguments

- moo:

  multiOmicDataSet object

- ...:

  additional arguments to pass to edgeR::cpm()

## Value

multiOmicDataSet with cpm-transformed counts

## See also

Other moo methods:
[`batch_correct_counts()`](https://ccbr.github.io/MOSuite/reference/batch_correct_counts.md),
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
[`set_color_pal()`](https://ccbr.github.io/MOSuite/reference/set_color_pal.md),
[`set_default_colors()`](https://ccbr.github.io/MOSuite/reference/set_default_colors.md)

## Examples

``` r
sample_meta <- data.frame(
  sample_id = c("KO_S3", "KO_S4", "WT_S1", "WT_S2"),
  condition = factor(
    c("knockout", "knockout", "wildtype", "wildtype"),
    levels = c("wildtype", "knockout")
  )
)
moo <- create_multiOmicDataSet_from_dataframes(sample_meta, gene_counts) |>
  calc_cpm()
head(moo@counts$cpm)
#>              gene_id KO_S3 KO_S4 WT_S1 WT_S2
#> 1 ENSG00000121410.11     0     0     0     0
#> 2  ENSG00000268895.5     0     0     0     0
#> 3 ENSG00000148584.15     0     0     0     0
#> 4 ENSG00000175899.14     0     0     0     0
#> 5  ENSG00000245105.3     0     0     0     0
#> 6 ENSG00000166535.20     0     0     0     0
```
