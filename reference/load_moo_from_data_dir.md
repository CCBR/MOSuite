# Load multiOmicDataSet from data directory

Searches recursively for a single .rds file and loads its
multiOmicDataSet. An error is raised if no files or multiple files are
found, or if the object is not a multiOmicDataSet.

## Usage

``` r
load_moo_from_data_dir(data_dir = file.path("..", "data"))
```

## Arguments

- data_dir:

  path to data directory containing .rds file (default: `../data`)

## Value

loaded multiOmicDataSet object

## Examples

``` r
if (FALSE) { # \dontrun{
moo <- load_moo_from_data_dir()
} # }
```
