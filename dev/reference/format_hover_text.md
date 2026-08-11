# Format hover text for interactive plots

Builds a character vector of hover labels from one or two metadata
columns. Validates that required columns exist in `plot_data` and raises
an error if any are missing.

## Usage

``` r
format_hover_text(
  plot_data,
  primary_colname,
  secondary_colname = NULL,
  missing_col_context = "plot",
  require_secondary = TRUE
)
```

## Arguments

- plot_data:

  data frame containing plot metadata columns.

- primary_colname:

  name of the primary column to display in hover text.

- secondary_colname:

  optional name of a secondary column appended below the primary label.
  If `NULL` or not present in `plot_data`, only the primary label is
  returned.

- missing_col_context:

  string describing the plot context, used in the error message when
  required columns are absent.

- require_secondary:

  if `TRUE` and `secondary_colname` is not `NULL`, the secondary column
  is also validated as required.

## Value

character vector of hover text labels, one per row of `plot_data`.
