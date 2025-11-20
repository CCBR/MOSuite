# Plot histogram for multiOmicDataSet

Plot histogram for multiOmicDataSet

## Arguments

- moo_counts:

  counts dataframe or `multiOmicDataSet` containing `count_type` &
  `sub_count_type` in the counts slot

- count_type:

  Required if `moo_counts` is a `multiOmicDataSet`: the type of counts
  to use – must be a name in the counts slot (`moo@counts`).

- sub_count_type:

  Used if `moo_counts` is a `multiOmicDataSet` AND if `count_type` is a
  list, specify the sub count type within the list

- ...:

  arguments forwarded to method:
  [plot_histogram_dat](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram_dat.md)

## See also

[plot_histogram](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram.md)
generic

Other plotters for multiOmicDataSets:
[`plot_corr_heatmap_moo`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap_moo.md),
[`plot_expr_heatmap_moo`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap_moo.md),
[`plot_pca_moo`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca_moo.md),
[`plot_read_depth_moo`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth_moo.md)

## Examples

``` r
# plot histogram for a counts slot in a multiOmicDataset Object
moo <- multiOmicDataSet(
  sample_metadata = nidap_sample_metadata,
  anno_dat = data.frame(),
  counts_lst = list("raw" = nidap_raw_counts)
)
p <- plot_histogram(moo, count_type = "raw")

# customize the plot
plot_histogram(moo,
  count_type = "raw",
  group_colname = "Group", color_by_group = TRUE
)

```
