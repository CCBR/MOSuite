# Compute colour legend text size

Computes legend text size from legend labels. Short legends keep the
larger default text used by simple group legends, while longer or denser
legends are scaled down.

## Usage

``` r
get_legend_text_size(
  labels,
  legend_text_size = NULL,
  min_legend_text_size = 8,
  max_legend_text_size = 18
)
```

## Arguments

- labels:

  Character vector of legend labels.

- legend_text_size:

  Optional explicit legend text size. When supplied, this value is
  returned unchanged.

- min_legend_text_size:

  Smallest automatically selected legend text size.

- max_legend_text_size:

  Largest automatically selected legend text size.

## Value

Numeric legend text size.
