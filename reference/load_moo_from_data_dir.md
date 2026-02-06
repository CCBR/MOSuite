# Load multiOmicDataSet from data directory

Searches the ../data directory for .rds files and loads the first
matching multiOmicDataSet object. Validates that the loaded object is of
the correct class.

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
