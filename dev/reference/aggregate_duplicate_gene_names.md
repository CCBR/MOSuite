# Aggregate duplicate gene names

Aggregate duplicate gene names

## Usage

``` r
aggregate_duplicate_gene_names(
  counts_dat,
  gene_name_column_to_use_for_collapsing_duplicates,
  aggregate_rows_with_duplicate_gene_names,
  split_gene_name
)
```

## Arguments

- counts_dat:

  dataframe with raw counts data

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

- split_gene_name:

  If `TRUE`, split the gene name column by any of these special
  characters: `,|_-:`

## Value

data frame with columns separated if possible
