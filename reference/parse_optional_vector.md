# Parse comma-separated string into a vector

Splits a comma-separated string into a trimmed character vector. Returns
NULL if input is empty, NULL, or has zero length.

## Usage

``` r
parse_optional_vector(x)
```

## Arguments

- x:

  character string with comma-separated values

## Value

character vector or NULL if input is empty

## Examples

``` r
if (FALSE) { # \dontrun{
parse_optional_vector("a, b, c")
parse_optional_vector("")
} # }
```
