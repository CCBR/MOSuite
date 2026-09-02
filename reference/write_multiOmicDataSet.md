# Write a multiOmicDataSet to RDS

Re-exported from
[MOObject::write_multiOmicDataSet](https://ccbr.github.io/MOObject/reference/write_multiOmicDataSet.html)
– view the MOObject docs for details.

## Usage

``` r
write_multiOmicDataSet(moo, filepath = "moo.rds")
```

## Arguments

- moo:

  [multiOmicDataSet](https://ccbr.github.io/MOObject/reference/multiOmicDataSet.html)
  object to serialize

- filepath:

  Path to the RDS file to write (default: "moo.rds")

## Value

Invisibly returns `filepath`.

## See also

[`MOObject::write_multiOmicDataSet()`](https://ccbr.github.io/MOObject/reference/write_multiOmicDataSet.html)

Other moo IO:
[`create_multiOmicDataSet_from_dataframes()`](https://ccbr.github.io/MOSuite/reference/create_multiOmicDataSet_from_dataframes.md),
[`create_multiOmicDataSet_from_files()`](https://ccbr.github.io/MOSuite/reference/create_multiOmicDataSet_from_files.md),
[`multiOmicDataSet()`](https://ccbr.github.io/MOSuite/reference/multiOmicDataSet.md),
[`read_multiOmicDataSet()`](https://ccbr.github.io/MOSuite/reference/read_multiOmicDataSet.md),
[`read_multiOmicDataSet_properties()`](https://ccbr.github.io/MOSuite/reference/read_multiOmicDataSet_properties.md),
[`write_multiOmicDataSet_properties()`](https://ccbr.github.io/MOSuite/reference/write_multiOmicDataSet_properties.md)
