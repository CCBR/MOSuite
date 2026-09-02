# Print CLI usage information

Writes a usage summary for the `mosuite` command-line tool to a
connection, typically
[`stderr()`](https://rdrr.io/r/base/showConnections.html).

## Usage

``` r
cli_usage(con = stderr())
```

## Arguments

- con:

  connection to write usage text to. Defaults to
  [`stderr()`](https://rdrr.io/r/base/showConnections.html).

## Value

invisibly returns `NULL`.
