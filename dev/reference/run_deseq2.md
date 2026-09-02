# Run DESeq2 on a multiOmicDataSet

Run DESeq2 on a multiOmicDataSet

## Usage

``` r
run_deseq2(moo, design, ...)
```

## Arguments

- moo:

  multiOmicDataSet object

- design:

  model formula for experimental design. Columns must exist in
  `meta_dat`.

- ...:

  remaining variables are forwarded to
  [`DESeq2::DESeq()`](https://rdrr.io/pkg/DESeq2/man/DESeq.html).

## Value

multiOmicDataSet object with DESeq2 slot filled

## See also

Other moo methods:
[`batch_correct_counts()`](https://ccbr.github.io/MOSuite/dev/reference/batch_correct_counts.md),
[`calc_cpm()`](https://ccbr.github.io/MOSuite/dev/reference/calc_cpm.md),
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/dev/reference/clean_raw_counts.md),
[`diff_counts()`](https://ccbr.github.io/MOSuite/dev/reference/diff_counts.md),
[`extract_counts()`](https://ccbr.github.io/MOSuite/dev/reference/extract_counts.md),
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
[`set_color_pal()`](https://ccbr.github.io/MOSuite/dev/reference/set_color_pal.md),
[`set_default_colors()`](https://ccbr.github.io/MOSuite/dev/reference/set_default_colors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
moo <- create_multiOmicDataSet_from_files(
  system.file("extdata", "sample_metadata.tsv.gz",
    package = "MOSuite"
  ),
  system.file("extdata",
    "RSEM.genes.expected_count.all_samples.txt.gz",
    package = "MOSuite"
  )
) |> filter_counts()
moo <- run_deseq2(moo, ~condition)
} # }
```
