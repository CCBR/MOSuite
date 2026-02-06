# Set up capsule environment and directories

Initializes the results directory structure and logs installed R package
versions. This is a common setup task used across all Code Ocean
capsules.

## Usage

``` r
setup_capsule_environment(base_results_dir = file.path("..", "results"))
```

## Arguments

- base_results_dir:

  base path to results directory (default: `../results`)

## Value

invisibly returns a list with `results_dir` and `plots_dir` paths

## Examples

``` r
if (FALSE) { # \dontrun{
setup_capsule_environment()
} # }
```
