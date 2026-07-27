# Display colors for a multiOmicDataSet object

Plots a palette strip for each group column stored in
`moo@analyses$colors`, stacked vertically. Each strip shows the assigned
hex colors and their codes.

## Usage

``` r
display_colors(moo)
```

## Arguments

- moo:

  A `multiOmicDataSet` object (see
  [`create_multiOmicDataSet_from_dataframes()`](https://ccbr.github.io/MOSuite/reference/create_multiOmicDataSet_from_dataframes.md)).

## Value

A
[patchwork](https://patchwork.data-imaginist.com/reference/wrap_plots.html)
of
[ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
objects, one per group column in `moo@analyses$colors`.

## Examples

``` r
moo <- create_multiOmicDataSet_from_dataframes(nidap_sample_metadata, nidap_raw_counts)
display_colors(moo)
```
