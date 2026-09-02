# Plot read depth for multiOmicDataSet

Plot read depth for multiOmicDataSet

## Arguments

- count_type:

  the type of counts to use. Must be a name in the counts slot
  (`names(moo@counts)`).

- sub_count_type:

  used if `count_type` is a list in the counts slot: specify the sub
  count type within the list. Must be a name in
  `names(moo@counts[[count_type]])`.

- sample_id_colname:

  column in sample metadata containing sample IDs.

- group_colname:

  sample metadata column used to color bars. Leave blank to use the
  current single-color bar fill.

- color_values:

  colors used when `group_colname` is supplied. Named vectors are
  matched to group values; unnamed vectors follow group order and are
  extended with MOSuite colors when too few colors are supplied.
  Defaults to `NULL`; when `NULL`, `mosuite_palette` is used for
  `data.frame` dispatch and stored colors are used for
  `multiOmicDataSet` dispatch.

## Value

ggplot barplot

## See also

[`plot_read_depth()`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.md)
generic

Other plotters for multiOmicDataSets:
[`plot_corr_heatmap,MOObject::multiOmicDataSet-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap-multiOmicDataSet.md),
[`plot_histogram,MOObject::multiOmicDataSet-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram.multiOmicDataSet.md),
[`plot_pca,MOObject::multiOmicDataSet-method`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.multiOmicDataSet.md)

## Examples

``` r
# multiOmicDataSet
moo <- multiOmicDataSet(
  sample_metadata = nidap_sample_metadata,
  anno_dat = data.frame(),
  counts_lst = list(
    "raw" = nidap_raw_counts,
    "clean" = nidap_clean_raw_counts
  )
)

plot_read_depth(moo, count_type = "clean")

```
