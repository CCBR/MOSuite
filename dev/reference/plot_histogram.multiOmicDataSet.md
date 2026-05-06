# Plot histogram for multiOmicDataSet

Plot histogram for multiOmicDataSet

## Arguments

- count_type:

  Required if `moo_counts` is a `multiOmicDataSet`: the type of counts
  to use – must be a name in the counts slot (`moo@counts`).

- sub_count_type:

  Used if `moo_counts` is a `multiOmicDataSet` AND if `count_type` is a
  list, specify the sub count type within the list

## See also

[`plot_histogram()`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram.md)
generic

Other plotters for multiOmicDataSets:
[`plot_corr_heatmap,MOSuite::multiOmicDataSet-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap-multiOmicDataSet.md),
[`plot_pca,MOSuite::multiOmicDataSet-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.multiOmicDataSet.md),
[`plot_read_depth,MOSuite::multiOmicDataSet-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.multiOmicDataSet.md)

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
