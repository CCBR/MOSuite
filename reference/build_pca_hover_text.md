# Build hover text for PCA interactive plots

Constructs a character vector of hover labels for 2D and 3D PCA plots.
Delegates to
[`format_hover_text()`](https://ccbr.github.io/MOSuite/reference/format_hover_text.md)
using `label_colname` (or `sample_id_colname` when `label_colname` is
`NULL`) as the primary label and `group_colname` as the secondary label.

## Usage

``` r
build_pca_hover_text(
  pca_data,
  sample_id_colname,
  group_colname,
  label_colname = NULL
)
```

## Arguments

- pca_data:

  data frame of PCA results merged with sample metadata.

- sample_id_colname:

  column name containing sample identifiers.

- group_colname:

  column name containing sample group labels.

- label_colname:

  optional column name to use as the primary hover label instead of
  `sample_id_colname`. If `NULL`, `sample_id_colname` is used.

## Value

character vector of hover text, one entry per row of `pca_data`.
