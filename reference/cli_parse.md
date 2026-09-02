# Parse a single CLI argument string

Converts a raw CLI argument string to an R value. Logical-like strings
(`"true"`, `"false"`, etc.) are returned as `TRUE`/`FALSE`. Other
strings are parsed as R expressions; if the result is a language object,
the original string is returned as-is.

## Usage

``` r
cli_parse(text)
```

## Arguments

- text:

  character string to parse.

## Value

parsed R value, or `text` unchanged if parsing yields a language object.
