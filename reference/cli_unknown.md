# Build an unknown-command error message

Constructs an error message for an unrecognised CLI method name,
optionally suggesting similar exported function names based on edit
distance.

## Usage

``` r
cli_unknown(method, exports)
```

## Arguments

- method:

  the unrecognised method name supplied by the user.

- exports:

  character vector of exported function names from MOSuite.

## Value

character string with the error message.
