# Write multiOmicDataSet properties to disk as CSV files

Writes the properties of a multiOmicDataSet object to disk as separate
files in output_dir. Properties that are data frames are saved as CSV
files, while all other objects are saved as RDS files.

## Usage

``` r
write_multiOmicDataSet_properties(moo, output_dir = "moo")
```

## Arguments

- moo:

  `multiOmicDataSet` object to write properties from

- output_dir:

  Directory where the properties will be saved (default: "moo")

## Value

Invisibly returns the `output_dir` where the files were saved
