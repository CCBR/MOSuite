# Call an MOSuite function with arguments specified in a json file

Call an MOSuite function with arguments specified in a json file

## Usage

``` r
cli_from_json(method, json, debug = FALSE)
```

## Arguments

- method:

  function in MOSuite to call

- json:

  path to a JSON file containing arguments for the function.
  Additionally, the JSON can contain the following keys:

  - `moo_input_rds` - filepath to an existing MultiOmicsDataset object
    in RDS format. This is required if the MOSuite function contains
    `moo` as an argument.

  - `moo_output_rds` - filepath to write the result to.

- debug:

  when TRUE, do not call the command, just return the expression.

## Value

invisible returns the function call
