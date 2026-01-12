# Clean Raw Counts

This function checks the input raw counts matrix for common formatting
problems with feature identifiers and sample names. If feature IDs
contain multiple IDs separated by special characters (\| - , or space)
they will be split into multiple columns. If duplicate feature IDs are
detected the counts are summed across duplicate feature ID rows within
each sample. Invalid sample names will also be reported and can be
automatically corrected. If your sample names are corrected here, be
sure to make equivalent changes to your metadata table.

## Usage

``` r
clean_raw_counts(
  moo,
  count_type = "raw",
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  samples_to_rename = "",
  cleanup_column_names = TRUE,
  split_gene_name = TRUE,
  aggregate_rows_with_duplicate_gene_names = TRUE,
  gene_name_column_to_use_for_collapsing_duplicates = "",
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = "clean"
)
```

## Arguments

- moo:

  multiOmicDataSet object (see
  [`create_multiOmicDataSet_from_dataframes()`](https://ccbr.github.io/MOSuite/reference/create_multiOmicDataSet_from_dataframes.md))

- count_type:

  the type of counts to use – must be a name in the counts slot
  (`moo@counts`)

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

- samples_to_rename:

  If you do not have a Plot Labels Column in your sample metadata table,
  you can use this parameter to rename samples manually for display on
  the PCA plot. Use "Add item" to add each additional sample for
  renaming. Use the following format to describe which old name (in your
  sample metadata table) you want to rename to which new name: old_name:
  new_name

- cleanup_column_names:

  Invalid raw counts column names can cause errors in the downstream
  analysis. If this is `TRUE`, any invalid column names will be
  automatically altered to a correct format. These format changes will
  include adding an "X" as the first character in any column name that
  began with a numeral and replacing some special characters ("-,:. ")
  with underscores ("\_"). Invalid sample names and any changes made
  will be detailed.

- split_gene_name:

  If `TRUE`, split the gene name column by any of these special
  characters: `,|_-:`

- aggregate_rows_with_duplicate_gene_names:

  If a Feature ID (from the "Cleanup Column Names" parameter above) is
  found to be duplicated on multiple rows of the raw counts, the Log
  will report these Feature IDs. Using the default behavior (`TRUE`),
  the counts for all rows with a duplicate Feature IDs are aggregated
  into a single row. Counts are summed across duplicate Feature ID rows
  within each sample. Additional identifier columns, if present (e.g.
  Ensembl IDs), will be preserved and multiple matching identifiers in
  such additional columns will appear as comma-separated values in an
  aggregated row.

- gene_name_column_to_use_for_collapsing_duplicates:

  Select the column with Feature IDs to use as grouping elements to
  collapse the counts matrix. The log output will list the columns
  available to identify duplicate row IDs in order to aggregate
  information. If left blank your "Feature ID" Column will be used to
  Aggregate Rows. If "Feature ID" column can be split into multiple IDs
  the non Ensembl ID name will be used to aggregate duplicate IDs. If
  "Feature ID" column does not contain Ensembl IDs the split Feature IDs
  will be named 'Feature_id_1' and 'Feature_id_2'. For this case an
  error will occur and you will have to manually enter the Column ID for
  this field.

- print_plots:

  Whether to print plots during analysis (Defaults to `FALSE`,
  overwritable using option 'moo_print_plots' or environment variable
  'MOO_PRINT_PLOTS')

- save_plots:

  Whether to save plots to files during analysis (Defaults to `TRUE`,
  overwritable using option 'moo_save_plots' or environment variable
  'MOO_SAVE_PLOTS')

- plots_subdir:

  subdirectory in `figures/` where plots will be saved if `save_plots`
  is `TRUE`

## Value

`multiOmicDataSet` with cleaned counts

## See also

Other moo methods:
[`batch_correct_counts()`](https://ccbr.github.io/MOSuite/reference/batch_correct_counts.md),
[`diff_counts()`](https://ccbr.github.io/MOSuite/reference/diff_counts.md),
[`filter_counts()`](https://ccbr.github.io/MOSuite/reference/filter_counts.md),
[`filter_diff()`](https://ccbr.github.io/MOSuite/reference/filter_diff.md),
[`normalize_counts()`](https://ccbr.github.io/MOSuite/reference/normalize_counts.md),
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_corr_heatmap.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/reference/plot_expr_heatmap.md),
[`plot_histogram()`](https://ccbr.github.io/MOSuite/reference/plot_histogram.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/reference/plot_pca.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/reference/plot_read_depth.md),
[`run_deseq2()`](https://ccbr.github.io/MOSuite/reference/run_deseq2.md),
[`set_color_pal()`](https://ccbr.github.io/MOSuite/reference/set_color_pal.md)

## Examples

``` r
moo <- create_multiOmicDataSet_from_dataframes(
  as.data.frame(nidap_sample_metadata),
  as.data.frame(nidap_raw_counts),
  sample_id_colname = "Sample",
) %>%
  clean_raw_counts(sample_id_colname = "Sample", feature_id_colname = "GeneName")
#> Saving 6.67 x 6.67 in image
#> * cleaning raw counts
#> Not able to identify multiple id's in GeneName
#> Columns that can be used to aggregate gene information GeneName
#> Aggregating the counts for the same ID in different chromosome locations.
#> Column used to Aggregate duplicate IDs: GeneName
#> Number of rows before Collapse: 43280
#> no duplicated IDs in GeneName
head(moo@counts$clean)
#>        GeneName A1 A2 A3 B1 B2 B3 C1 C2 C3
#> 1 RP23-271O17.1  0  0  0  0  0  0  0  0  0
#> 2       Gm26206  0  0  0  0  0  0  0  0  0
#> 3          Xkr4  0  0  0  0  0  0  0  0  0
#> 4 RP23-317L18.1  0  0  0  0  0  0  0  0  0
#> 5 RP23-317L18.4  0  0  0  0  0  0  0  0  0
#> 6 RP23-317L18.3  0  0  0  0  0  0  0  0  0
```
