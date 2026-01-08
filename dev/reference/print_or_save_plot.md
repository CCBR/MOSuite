# Print and/or save a ggplot

If `save_plots` is `TRUE`, the plot will be saved as an image to the
path at `file.path(plots_dir, filename)`. If `plot_obj` is a ggplot,
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
is used to save the image. Otherwise, `graphics_device` is used
(`grDevice::png()` by default).

## Usage

``` r
print_or_save_plot(
  plot_obj,
  filename,
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_dir = options::opt("plots_dir"),
  graphics_device = grDevices::png,
  ...
)
```

## Arguments

- plot_obj:

  plot object (e.g. ggplot, ComplexHeatmap...)

- filename:

  name of the output file. will be joined with the `plots_dir` option.

- print_plots:

  Whether to print plots during analysis (Defaults to `FALSE`,
  overwritable using option 'moo_print_plots' or environment variable
  'MOO_PRINT_PLOTS')

- save_plots:

  Whether to save plots to files during analysis (Defaults to `TRUE`,
  overwritable using option 'moo_save_plots' or environment variable
  'MOO_SAVE_PLOTS')

- plots_dir:

  Path where plots are saved when `moo_save_plots` is `TRUE` (Defaults
  to `"figures/"`, overwritable using option 'moo_plots_dir' or
  environment variable 'MOO_PLOTS_DIR')

- graphics_device:

  Default: `grDevice::png()`. Only used if the plot is not a ggplot.

- ...:

  arguments forwarded to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)

## Value

invisibly returns the path where the plot image was saved to the disk

## See also

Other plotters:
[`plot_corr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_corr_heatmap.md),
[`plot_expr_heatmap()`](https://ccbr.github.io/MOSuite/dev/reference/plot_expr_heatmap.md),
[`plot_histogram()`](https://ccbr.github.io/MOSuite/dev/reference/plot_histogram.md),
[`plot_pca()`](https://ccbr.github.io/MOSuite/dev/reference/plot_pca.md),
[`plot_read_depth()`](https://ccbr.github.io/MOSuite/dev/reference/plot_read_depth.md)
