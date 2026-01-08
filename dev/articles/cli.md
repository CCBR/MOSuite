# Calling MOSuite from the CLI

> ⚠️ **Most users do not need to use the CLI.** We recommend using
> MOSuite within R scripts, R Markdown, or Quarto documents for the vast
> majority of use-cases, as shown in the [**introductory
> vignette**](https://ccbr.github.io/MOSuite/articles/intro.html). The
> CLI is provided for a very specialized situation where MOSuite is run
> in an environment that cannot use R scripts natively.

MOSuite includes an executable file called `mosuite`. Any user-facing
function in the MOSuite R package can be called with
`mosuite [function]` from the unix CLI. Function arguments are passed in
via a JSON file. In addition to arguments used by the function, the JSON
file can contain the following keys:

- `moo_input_rds` - file path to an existing MultiOmicsDataset object in
  RDS format. This is required if the MOSuite function has `moo` as an
  argument (most user-facing functions do).
- `moo_output_rds` - file path to write the result to.

## Usage

Run `mosuite --help` in a unix shell to see the full CLI usage:

``` sh
Usage: mosuite [function] [--json=path/to/args.json]

[function] should be the name of a function exported from MOSuite.
[--json] should specify the path to a JSON file with arguments accepted by that function.
         The equals sign (=) is required to separate --json from the path.

Additionally, the JSON file can contain the following keys:
  - moo_input_rds: file path to an existing MultiOmicsDataset object in RDS format.
    This is required if `method` has `moo` as an argument.
  - moo_output_rds: file path to write the result to.

Use `mosuite [function] --help` for more information about the associated function.

Main functions:
  mosuite create_multiOmicDataSet_from_files
  mosuite filter_counts
  mosuite clean_raw_counts
  mosuite normalize_counts
  mosuite batch_correct_counts
  mosuite diff_counts
  mosuite filter_diff
```

## Installing the MOSuite CLI

### Docker Container

We provide a docker container with the MOSuite R package and CLI
installed as of v0.2.0 and later.
<https://hub.docker.com/r/nciccbr/mosuite>

Running this container with docker or singularity is the recommend way
to run MOSuite in pipelines and HPC environments.

``` sh
singularity exec docker://nciccbr/mosuite:v0.2.0 bash mosuite --help
singularity exec docker://nciccbr/mosuite:v0.2.0 R -s -e \
  'cat("MOSuite version:", installed.packages()["MOSuite",][["Version"]])'
```

### Installation on a personal computer

After installing the R package, you can use
[`system.file()`](https://rdrr.io/r/base/system.file.html) to locate the
`mosuite` executable file with R:

``` r
# remotes::install_github("CCBR/MOSuite", dependencies = TRUE)
system.file("exec", "mosuite", package = "MOSuite")
#> [1] "/home/runner/work/_temp/Library/MOSuite/exec/mosuite"
```

You should add this executable to your `PATH` environment variable.

``` sh
export PATH="$PATH:/path/to/exec/mosuite"
```

If you’re using the [MOSuite docker container](#docker-container), it is
already included in the path.

## Example end-to-end script

You can create a shell script to run the full MOSuite pipeline. This
script assumes you have a directory `json_args/` with JSON files to set
each function’s arguments.

``` bash
#!/usr/bin/env bash
set -euo pipefail

# set MOSuite options for plots
export MOO_SAVE_PLOTS=TRUE
export MOO_PLOTS_DIR=./figures
mkdir -p $MOO_PLOTS_DIR

# add mosuite executable to the path
mosuite=$(R -s -e "cat(system.file('exec','mosuite', package='MOSuite'))")
export PATH="$PATH:$(dirname $mosuite)"

mosuite create_multiOmicDataSet_from_files --json=json_args/common/create_multiOmicDataSet_from_files.json
mosuite clean_raw_counts --json=json_args/common/clean_raw_counts.json
mosuite filter_counts --json=json_args/common/filter_counts.json
mosuite normalize_counts --json=json_args/common/normalize_counts.json
mosuite batch_correct_counts --json=json_args/common/batch_correct_counts.json
mosuite diff_counts --json=json_args/common/diff_counts.json
mosuite filter_diff --json=json_args/common/filter_diff.json
```

The example script and accompanying JSON files are included in the
package data. You can copy them to your working directory with R:

``` r
# copy the example script
file.copy(
  system.file("extdata", "example_script.sh", package = "MOSuite"),
  to = "./"
)
# copy the JSON files
file.copy(
  system.file("extdata", "json_args", package = "MOSuite"),
  to = "./",
  recursive = TRUE
)
# copy the raw counts & sample metadata
file.copy(
  system.file("extdata", "nidap", "Raw_Counts.csv.gz", package = "MOSuite"),
  to = "./"
)
file.copy(
  system.file(
    "extdata",
    "nidap",
    "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz",
    package = "MOSuite"
  ),
  to = "./"
)
```

Then run the script from the CLI:

``` bash
bash ./example_script.sh
```

The final multiOmicDataSet will be in `moo.rds` and figures from each
step will be in `./figures/`.

## Writing JSON files

Create a JSON file with arguments for
[`create_multiOmicDataSet_from_files()`](https://ccbr.github.io/MOSuite/dev/reference/create_multiOmicDataSet_from_files.md).
You can use R code as below or write it by hand.

``` r
j <- list(
  feature_counts_filepath = system.file(
    "extdata",
    "RSEM.genes.expected_count.all_samples.txt.gz",
    package = "MOSuite"
  ),
  sample_meta_filepath = system.file(
    "extdata",
    "sample_metadata.tsv.gz",
    package = "MOSuite"
  ),
  moo_output_rds = "moo.rds"
)
jsonlite::write_json(j, "args_1.json")
```

In a unix shell, call
[`create_multiOmicDataSet_from_files()`](https://ccbr.github.io/MOSuite/dev/reference/create_multiOmicDataSet_from_files.md)
and specify the path to the JSON file:

``` bash
mosuite create_multiOmicDataSet_from_files --json=args_1.json
```

This is equivalent to running the following R code:

``` r
library(MOSuite)
moo <- create_multiOmicDataSet_from_files(
  feature_counts_filepath = system.file(
    "extdata",
    "RSEM.genes.expected_count.all_samples.txt.gz",
    package = "MOSuite"
  ),
  sample_meta_filepath = system.file(
    "extdata",
    "sample_metadata.tsv.gz",
    package = "MOSuite"
  )
)
#> Rows: 58929 Columns: 6
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (2): gene_id, GeneName
#> dbl (4): KO_S3, KO_S4, WT_S1, WT_S2
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 4 Columns: 2
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (2): sample_id, condition
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
readr::write_rds(moo, "moo.rds")
```

You can use the `moo` object you just created as input to other MOSuite
functions.

Create a JSON file of arguments for
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/dev/reference/clean_raw_counts.md)
with R (or write it by hand):

``` r
j <- list(
  moo_input_rds = "moo.rds",
  moo_output_rds = "moo.rds",
  save_plots = TRUE
)
jsonlite::write_json(j, "args_2.json")
```

Then run
[`clean_raw_counts()`](https://ccbr.github.io/MOSuite/dev/reference/clean_raw_counts.md):

``` bash
mosuite clean_raw_counts --json=args_2.json
```

Results are saved to `moo.rds`. Overwriting the same `moo` file is
recommended to save disk space, as the multiOmicDataset object saves
intermediate results within its data structure.

## Template JSON files

JSON file templates with default arguments for the main functions are
bundled with the package. You can copy them to your current directory
like so:

``` r
file.copy(
  system.file("extdata", "json_args", "defaults", package = "MOSuite"),
  to = "./",
  recursive = TRUE
)
#> [1] TRUE
```
