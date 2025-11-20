# MOSuite Options

Internally used, package-specific options. All options will prioritize R
options() values, and fall back to environment variables if undefined.
If neither the option nor the environment variable is set, a default
value is used.

## Checking Option Values

Option values specific to `MOSuite` can be accessed by passing the
package name to `env`.

    options::opts(env = "MOSuite")

    options::opt(x, default, env = "MOSuite")

## Options

- print_plots:

  default:

  :   FALSE

  option:

  :   moo_print_plots

  envvar:

  :   MOO_PRINT_PLOTS (evaluated if possible, raw string otherwise)

- save_plots:

  default:

  :   FALSE

  option:

  :   moo_save_plots

  envvar:

  :   MOO_SAVE_PLOTS (evaluated if possible, raw string otherwise)

- plots_dir:

  default:

  :   "figures/"

  option:

  :   moo_plots_dir

  envvar:

  :   MOO_PLOTS_DIR (evaluated if possible, raw string otherwise)

## See also

options getOption Sys.setenv Sys.getenv
