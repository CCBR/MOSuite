# Get constructor-style default colors for a multiOmicDataSet column

Returns stored colors from `moo@analyses$colors[[colname]]` when
available; otherwise recreates the default colors that MOSuite
constructors assigned before the object class moved to MOObject.

## Usage

``` r
get_moo_default_colors(moo, colname, palette = mosuite_palette)
```

## Arguments

- moo:

  `multiOmicDataSet` object

- colname:

  column name in `moo@sample_meta`

- palette:

  Character vector of colors to assign. Defaults to `mosuite_palette`.

## Value

Named character vector of colors, or `NULL` when `colname` is not
present in `moo@sample_meta`.
