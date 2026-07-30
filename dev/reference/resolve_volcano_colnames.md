# Resolve volcano plot column names

Auto-detects `change_colname` and `signif_colname` from a data frame
when either is `NULL`. Used by
[`plot_volcano_enhanced()`](https://ccbr.github.io/MOSuite/dev/reference/plot_volcano_enhanced.md)
and
[`plot_volcano_summary()`](https://ccbr.github.io/MOSuite/dev/reference/plot_volcano_summary.md).

## Usage

``` r
resolve_volcano_colnames(diff_dat, change_colname, signif_colname)
```

## Arguments

- diff_dat:

  A data frame of differential analysis results.

- change_colname:

  Character vector of logFC column names, or `NULL` to auto-detect
  columns ending in `_logFC`.

- signif_colname:

  Character vector of significance column names, or `NULL` to
  auto-detect, preferring `_adjpval` over `_pval`.

## Value

A named list with elements `change_colname` and `signif_colname`.

## Details

When `change_colname` is `NULL`, all columns ending in `_logFC` are
used. When `signif_colname` is `NULL`, significance columns are detected
by checking for `_adjpval` columns first, then `_pval`, for each
contrast derived from `change_colname`.
