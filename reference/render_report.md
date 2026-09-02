# Render the template report

Copy the Quarto template to the current working directory and render it
using
[`quarto::quarto_render()`](https://quarto-dev.github.io/quarto-r/reference/quarto_render.html).
The rendered report will be saved in the current working directory. You
can specify additional arguments to
[`quarto::quarto_render()`](https://quarto-dev.github.io/quarto-r/reference/quarto_render.html)
to customize the rendering process.

## Usage

``` r
render_report(
  qmd_template = system.file("quarto", "report.qmd", package = "MOSuite"),
  qmd_src = NULL,
  ...
)
```

## Arguments

- qmd_template:

  Path to the Quarto report file (default is the template report in the
  package).

- qmd_src:

  Optional path to copy the Quarto report template to before rendering.
  If `NULL` (default), the template will be copied to the current
  working directory with the same filename as the template. If a file
  already exists at `qmd_src`, it will not be overwritten.

- ...:

  Additional arguments passed to
  [`quarto::quarto_render()`](https://quarto-dev.github.io/quarto-r/reference/quarto_render.html),
  such as `execute_params` (a named list of parameters) or `quarto_args`
  (a character vector of CLI flags, e.g.
  `c("--output-dir", "/path/to/out")`).

## Details

You can edit the copy of `report.qmd` in the current working directory
to customize the report.

## Examples

``` r
if (FALSE) { # \dontrun{
render_report(execute_params = list(
  counts_csv = system.file("extdata", "nidap", "Raw_Counts.csv.gz",
                           package = "MOSuite"),
  samplesheet_csv = system.file("extdata", "nidap",
    "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz",
    package = "MOSuite")
))

# Render to a specific output directory
render_report(
  quarto_args = c("--output-dir", "./results"),
  execute_params = list(
    counts_csv = system.file("extdata", "nidap", "Raw_Counts.csv.gz",
                             package = "MOSuite"),
    samplesheet_csv = system.file("extdata", "nidap",
      "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz",
      package = "MOSuite")
  )
)
} # }
```
