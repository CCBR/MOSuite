# Construct a multiOmicDataSet object from text files (e.g. TSV, CSV).

Wraps
[`MOObject::create_multiOmicDataSet_from_files()`](https://ccbr.github.io/MOObject/reference/create_multiOmicDataSet_from_files.html)
to add a "colors" analysis if it is not already present.

## Usage

``` r
create_multiOmicDataSet_from_files(
  sample_meta_filepath,
  feature_counts_filepath,
  count_type = "raw",
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  delim = NULL,
  ...
)
```

## Arguments

- sample_meta_filepath:

  path to text file with sample IDs and metadata for differential
  analysis.

- feature_counts_filepath:

  path to text file of expected feature counts (e.g. gene counts from
  RSEM).

- count_type:

  type to assign the values of `counts_dat` to in the `counts` slot

- sample_id_colname:

  name of the column in `sample_metadata` that contains the sample IDs.
  (Default: `NULL` - first column in the sample metadata will be used.)

- feature_id_colname:

  name of the column in `counts_dat` that contains feature/gene IDs.
  (Default: `NULL` - first column in the count data will be used.)

- delim:

  Delimiter used in the input files. Any delimiter accepted by
  [`readr::read_delim()`](https://readr.tidyverse.org/reference/read_delim.html)
  can be used. If the files are in CSV format, set `delim = ','`; for
  TSV format, set `delim = '\t'`.

- ...:

  additional arguments forwarded to
  [`readr::read_delim()`](https://readr.tidyverse.org/reference/read_delim.html).

## Value

A
[MOObject::multiOmicDataSet](https://ccbr.github.io/MOObject/reference/multiOmicDataSet.html)
object.

## See also

[`MOObject::create_multiOmicDataSet_from_files()`](https://ccbr.github.io/MOObject/reference/create_multiOmicDataSet_from_files.html)

Other moo IO:
[`create_multiOmicDataSet_from_dataframes()`](https://ccbr.github.io/MOSuite/reference/create_multiOmicDataSet_from_dataframes.md),
[`multiOmicDataSet()`](https://ccbr.github.io/MOSuite/reference/multiOmicDataSet.md),
[`read_multiOmicDataSet()`](https://ccbr.github.io/MOSuite/reference/read_multiOmicDataSet.md),
[`read_multiOmicDataSet_properties()`](https://ccbr.github.io/MOSuite/reference/read_multiOmicDataSet_properties.md),
[`write_multiOmicDataSet()`](https://ccbr.github.io/MOSuite/reference/write_multiOmicDataSet.md),
[`write_multiOmicDataSet_properties()`](https://ccbr.github.io/MOSuite/reference/write_multiOmicDataSet_properties.md)
