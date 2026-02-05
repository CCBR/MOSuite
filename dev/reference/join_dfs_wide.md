# Join dataframes in named list to wide dataframe

The first column is assumed to be shared by all dataframes

## Usage

``` r
join_dfs_wide(df_list, join_fn = dplyr::left_join)
```

## Arguments

- df_list:

  named list of dataframes

- join_fn:

  join function to use (Default:
  [`dplyr::left_join`](https://dplyr.tidyverse.org/reference/mutate-joins.html))

## Value

wide dataframe

## Examples

``` r
dfs <- list(
  "a_vs_b" = data.frame(id = c("a1", "b2", "c3"), score = runif(3)),
  "b_vs_c" = data.frame(id = c("a1", "b2", "c3"), score = rnorm(3))
)
dfs |> join_dfs_wide()
#> Joining with `by = join_by(id)`
#>   id a_vs_b_score b_vs_c_score
#> 1 a1    0.5315735   -0.8267890
#> 2 b2    0.4936370   -1.5123997
#> 3 c3    0.7793086    0.9353632
```
