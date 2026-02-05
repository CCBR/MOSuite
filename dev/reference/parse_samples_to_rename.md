# Parse sample rename pairs from string

Parses a string containing sample rename pairs in format
"old:new,old2:new2" and returns a named list where names are old sample
names and values are new names.

## Usage

``` r
parse_samples_to_rename(x)
```

## Arguments

- x:

  character string with rename pairs in format "old:new,old2:new2"

## Value

named list with old names as keys and new names as values, or NULL if
empty

## Examples

``` r
if (FALSE) { # \dontrun{
parse_samples_to_rename("sample1:S1,sample2:S2")
parse_samples_to_rename("")
} # }
```
