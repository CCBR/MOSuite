# Print help for a CLI method

Displays the R help page for the named MOSuite function in the `MOSuite`
package.

## Usage

``` r
cli_help(method)
```

## Arguments

- method:

  name of the MOSuite function to show help for.

## Value

result of `print(utils::help(...))`, invisibly.
