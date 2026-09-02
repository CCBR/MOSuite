# Compute a wrapped colour legend column count

Computes a conservative number of legend columns for horizontal ggplot
colour legends. Top and bottom legends are wrapped based on the number
of labels and the longest label length. Other legend positions return
`NULL` so their existing ggplot layout is preserved.

## Usage

``` r
get_legend_column_count(
  labels,
  legend_position = "top",
  ncol = NULL,
  legend_text_size = 10,
  max_label_characters_per_row = 45
)
```

## Arguments

- labels:

  Character vector of legend labels.

- legend_position:

  Legend position passed to
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html).

- ncol:

  Optional maximum number of legend columns.

- legend_text_size:

  Legend text size used to scale the horizontal space estimate. Larger
  legend text uses fewer columns.

- max_label_characters_per_row:

  Approximate total label characters to fit on one horizontal legend
  row.

## Value

Integer column count for top/bottom legends, or `NULL` when no wrapping
should be applied.
