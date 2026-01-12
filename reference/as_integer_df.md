# Convert all numeric columns in a dataframe to integers

Round doubles to integers and convert to integer type

## Usage

``` r
as_integer_df(counts_tbl)
```

## Arguments

- counts_tbl:

  data frame with numeric columns

## Value

data frame with any numeric columns as integers

## Examples

``` r
if (FALSE) { # \dontrun{
data.frame(a = c(0, 0.1, 2.3, 5L, 6.9)) %>% as_integer_df()
} # }
```
