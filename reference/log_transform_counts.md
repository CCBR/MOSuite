# Log-transform count columns in a data frame

Log-transform count columns in a data frame

## Usage

``` r
log_transform_counts(
  counts_dat,
  feature_id_colname = NULL,
  sample_colnames = NULL,
  pseudocount = 0.5,
  base = "ln"
)
```

## Arguments

- counts_dat:

  data frame of feature counts.

- feature_id_colname:

  name of the column in `counts_dat` that contains feature/gene IDs.
  (Default: `NULL` - first column in the count data will be used.)

- sample_colnames:

  optional vector of sample columns to transform. If `NULL`, all columns
  except `feature_id_colname` are transformed.

- pseudocount:

  value added before log transformation.

- base:

  logarithm base to use for the transformation. Use a numeric value, or
  `"e"`, `"ln"`, or `"natural"` for natural log. Default is `"ln"`.

## Value

count data frame with selected count columns transformed as
`log(x + pseudocount, base)`.
