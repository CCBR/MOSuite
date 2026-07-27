# Remove low-count genes

TODO this function also transforms raw counts to CPM, but that should be
a separate function before this step, before filter_counts function()

## Usage

``` r
remove_low_count_genes(
  counts_dat,
  sample_metadata,
  sample_id_colname = NULL,
  feature_id_colname,
  group_colname,
  use_cpm_counts_to_filter = TRUE,
  use_group_based_filtering = FALSE,
  minimum_count_value_to_be_considered_nonzero = 8,
  minimum_number_of_samples_with_nonzero_counts_in_total = 7,
  minimum_number_of_samples_with_nonzero_counts_in_a_group = 3
)
```

## Arguments

- sample_id_colname:

  The column from the sample metadata containing the sample names. The
  names in this column must exactly match the names used as the sample
  column names of your input Counts Matrix. (Default: `NULL` - first
  column in the sample metadata will be used.)

- feature_id_colname:

  The column from the counts data containing the Feature IDs (Usually
  Gene or Protein ID). This is usually the first column of your input
  Counts Matrix. Only columns of Text type from your input Counts Matrix
  will be available to select for this parameter. (Default: `NULL` -
  first column in the counts matrix will be used.)

- group_colname:

  The column from the sample metadata containing the sample group
  information. This is usually a column showing to which experimental
  treatments each sample belongs (e.g. WildType, Knockout, Tumor,
  Normal, Before, After, etc.).

- use_cpm_counts_to_filter:

  If no transformation has been been performed on counts matrix (eg Raw
  Counts) set to TRUE. If TRUE counts will be transformed to CPM and
  filtered based on given criteria. If gene counts matrix has been
  transformed (eg log2, CPM, FPKM or some form of Normalization) set to
  FALSE. If FALSE no further transformation will be applied and features
  will be filtered as is. For RNAseq data RAW counts should be
  transformed to CPM in order to properly filter.

- use_group_based_filtering:

  If TRUE, only keeps features (e.g. genes) that have at least a certain
  number of samples passing the threshold in at least one group

- minimum_count_value_to_be_considered_nonzero:

  Minimum value in the selected filtering table required for a sample to
  be considered nonzero. If `use_cpm_counts_to_filter` is `TRUE`, this
  threshold is applied to CPM values. If `use_cpm_counts_to_filter` is
  `FALSE`, this threshold is applied directly to the selected
  `count_type` table.

- minimum_number_of_samples_with_nonzero_counts_in_total:

  Minimum number of samples in total that must meet the
  `minimum_count_value_to_be_considered_nonzero` threshold for a feature
  to be kept.

- minimum_number_of_samples_with_nonzero_counts_in_a_group:

  Only keeps genes that have at least this number of samples meeting the
  threshold in at least one group

## Value

counts matrix with low-count genes removed
