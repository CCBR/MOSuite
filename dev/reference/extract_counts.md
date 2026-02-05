# Extract count data

Extract count data

## Usage

``` r
extract_counts(moo, count_type, sub_count_type = NULL)
```

## Arguments

- moo:

  multiOmicDataSet containing `count_type` & `sub_count_type` in the
  counts slot

- count_type:

  the type of counts to use – must be a name in the counts slot
  (`moo@counts[[count_type]]`)

- sub_count_type:

  if `count_type` is a list, specify the sub count type within the list
  (`moo@counts[[count_type]][[sub_count_type]]`). (Default: `NULL`)

## Examples

``` r
moo <- multiOmicDataSet(
  sample_metadata = as.data.frame(nidap_sample_metadata),
  anno_dat = data.frame(),
  counts_lst = list(
    "raw" = as.data.frame(nidap_raw_counts),
    "clean" = as.data.frame(nidap_clean_raw_counts),
    "filt" = as.data.frame(nidap_filtered_counts),
    "norm" = list(
      "voom" = as.data.frame(nidap_norm_counts)
    )
  )
)

moo |>
  extract_counts("filt") |>
  head()
#>            Gene   A1  A2   A3   B1   B2  B3   C1  C2   C3
#> 1 0610007P14Rik 1049 950  934 1068 1140 947 1393 907 1427
#> 2 0610009B22Rik  283 590  615  241  383 608  299 186  696
#> 3 0610010F05Rik  352 678 1377  958  879 616  332   0  186
#> 4 0610011F06Rik  430 565  553  462  558 688  710 826  706
#> 5 0610012G03Rik  480 589  683  324  596 673  909 933  419
#> 6 0610037L13Rik  467 570  593  558  330 423  356 198  568

moo |>
  extract_counts("norm", "voom") |>
  head()
#>            Gene       A1       A2       A3       B1       B2       B3       C1
#> 1 0610007P14Rik 6.532994 6.192871 5.954869 6.375896 6.275880 6.119449 6.419913
#> 2 0610009B22Rik 4.484983 5.448875 5.286875 3.445612 4.451347 5.473886 3.500359
#> 3 0610010F05Rik 4.883688 5.668494 6.537590 6.216408 5.893089 5.498884 3.845207
#> 4 0610011F06Rik 5.199684 5.374085 5.112952 5.155558 5.163359 5.650929 5.441965
#> 5 0610012G03Rik 5.368118 5.445918 5.456511 4.567138 5.274928 5.625039 5.787457
#> 6 0610037L13Rik 5.327987 5.388747 5.233520 5.450169 3.656585 4.929386 4.274944
#>         C2       C3
#> 1 6.172204 6.497050
#> 2 4.709254 5.471951
#> 3 2.685177 2.805426
#> 4 6.043492 5.490958
#> 5 6.214163 4.682896
#> 6 4.744405 5.173531
```
