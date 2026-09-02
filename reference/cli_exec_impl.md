# Internal implementation of CLI execution

Parses command-line arguments and dispatches to the appropriate MOSuite
function. Called by
[`cli_exec()`](https://ccbr.github.io/MOSuite/reference/cli_exec.md).

## Usage

``` r
cli_exec_impl(clargs)
```

## Arguments

- clargs:

  character vector of command-line arguments.

## Value

result of the dispatched MOSuite function call.
