# Write a multiOmicDataSet to disk as an RDS file

Write a multiOmicDataSet to disk as an RDS file

## Usage

``` r
write_multiOmicDataSet(moo, filepath = "moo.rds")
```

## Arguments

- moo:

  [multiOmicDataSet](https://ccbr.github.io/MOSuite/reference/multiOmicDataSet.md)
  object to serialize

- filepath:

  Path to the RDS file to write (default: "moo.rds")

## Value

Invisibly returns `filepath`
