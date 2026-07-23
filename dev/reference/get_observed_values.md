# Get observed values from a column

Returns non-missing values from `dat[[colname]]`. For factor columns,
values are returned in factor-level order; otherwise, values keep
first-observed order.

## Usage

``` r
get_observed_values(dat, colname)
```

## Arguments

- dat:

  data frame

- colname:

  column name in `dat`

## Value

character vector of observed values
