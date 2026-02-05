# Bind dataframes in named list to long dataframe

The dataframes must have all of the same columns

## Usage

``` r
bind_dfs_long(df_list, outcolname = contrast)
```

## Arguments

- df_list:

  named list of dataframes

- outcolname:

  column name in output dataframe for the names from the named list

## Value

long dataframe with new column `outcolname` from named list

## Examples

``` r
dfs <- list(
  "a_vs_b" = data.frame(id = c("a1", "b2", "c3"), score = runif(3)),
  "b_vs_c" = data.frame(id = c("a1", "b2", "c3"), score = rnorm(3))
)
dfs |> bind_dfs_long()
#>   id contrast       score
#> 1 a1   a_vs_b  0.87460066
#> 2 b2   a_vs_b  0.17494063
#> 3 c3   a_vs_b  0.03424133
#> 4 a1   b_vs_c -0.46662045
#> 5 b2   b_vs_c -0.85719040
#> 6 c3   b_vs_c -1.52474418
```
