# Parse comma-separated string with default fallback

Splits a comma-separated string into a trimmed character vector. Returns
a default value if input is empty, NULL, or has zero length.

## Usage

``` r
parse_vector_with_default(x, default)
```

## Arguments

- x:

  character string with comma-separated values

- default:

  default value to return if x is empty

## Value

character vector or default value

## Examples

``` r
if (FALSE) { # \dontrun{
parse_vector_with_default("a, b, c", "default")
parse_vector_with_default("", "default")
} # }
```
