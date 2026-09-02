# Add wrapped colour legend layout to a ggplot

Applies a colour guide with a wrapped column count for top and bottom
legends. Left, right, and hidden legends are returned unchanged.

## Usage

``` r
add_colour_legend_layout(
  plot,
  labels,
  legend_position = "top",
  ncol = NULL,
  legend_text_size = 10,
  max_label_characters_per_row = 45,
  guide_override_aes = NULL
)
```

## Arguments

- plot:

  A `ggplot2` plot object.

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

A `ggplot2` plot object with colour legend layout applied when needed.
