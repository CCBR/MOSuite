# Construct a multiOmicDataSet object from data frames

Wraps
[`MOObject::create_multiOmicDataSet_from_dataframes()`](https://ccbr.github.io/MOObject/reference/create_multiOmicDataSet_from_dataframes.html)
to add a "colors" analysis if it is not already present.

## Usage

``` r
create_multiOmicDataSet_from_dataframes(
  sample_metadata,
  counts_dat,
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  count_type = "raw"
)
```

## Arguments

- sample_metadata:

  sample metadata as a data frame or tibble. The first column is assumed
  to contain the sample IDs which must correspond to column names in the
  raw counts.

- counts_dat:

  data frame of feature counts (e.g. expected feature counts from RSEM).

- sample_id_colname:

  name of the column in `sample_metadata` that contains the sample IDs.
  (Default: `NULL` - first column in the sample metadata will be used.)

- feature_id_colname:

  name of the column in `counts_dat` that contains feature/gene IDs.
  (Default: `NULL` - first column in the count data will be used.)

- count_type:

  type to assign the values of `counts_dat` to in the `counts` slot

## Value

A
[MOObject::multiOmicDataSet](https://ccbr.github.io/MOObject/reference/multiOmicDataSet.html)
object.

## See also

[`MOObject::create_multiOmicDataSet_from_dataframes()`](https://ccbr.github.io/MOObject/reference/create_multiOmicDataSet_from_dataframes.html)

Other moo IO:
[`create_multiOmicDataSet_from_files()`](https://ccbr.github.io/MOSuite/reference/create_multiOmicDataSet_from_files.md),
[`multiOmicDataSet()`](https://ccbr.github.io/MOSuite/reference/multiOmicDataSet.md),
[`read_multiOmicDataSet()`](https://ccbr.github.io/MOSuite/reference/read_multiOmicDataSet.md),
[`read_multiOmicDataSet_properties()`](https://ccbr.github.io/MOSuite/reference/read_multiOmicDataSet_properties.md),
[`write_multiOmicDataSet()`](https://ccbr.github.io/MOSuite/reference/write_multiOmicDataSet.md),
[`write_multiOmicDataSet_properties()`](https://ccbr.github.io/MOSuite/reference/write_multiOmicDataSet_properties.md)
