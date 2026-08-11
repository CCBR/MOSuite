# Resolve the log transform base to a numeric value

Converts string aliases (`"e"`, `"ln"`, `"natural"`) to `exp(1)` and
validates that numeric bases are positive and not equal to 1.

## Usage

``` r
resolve_log_transform_base(base)
```

## Arguments

- base:

  a single numeric value or one of `"e"`, `"ln"`, `"natural"`.

## Value

numeric base value suitable for use in `log(x, base)`.
