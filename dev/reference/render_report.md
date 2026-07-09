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
  such as a named list of parameters.

## Details

You can edit the copy of `report.qmd` in the current working directory
to customize the report.

## Examples

``` r
render_report(execute_params = list(
  counts_csv = system.file("extdata", "nidap", "Raw_Counts.csv.gz",
                           package = "MOSuite"),
  samplesheet_csv = system.file("extdata", "nidap",
    "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz",
    package = "MOSuite")
))
#> 
#> 
#> processing file: report.qmd
#> 1/17                   
#> 2/17 [deps]            
#> 3/17                   
#> 4/17 [initialize]      
#> 5/17                   
#> 6/17 [analyze]         
#> 7/17                   
#> 8/17 [pca_3D]          
#> 9/17                   
#> 10/17 [expr_heatmap]    
#> 11/17                   
#> 12/17 [volcano_summary] 
#> 13/17                   
#> 14/17 [volcano_enhanced]
#> 15/17                   
#> 16/17 [venn_diagram]    
#> 17/17                   
#> output file: report.knit.md
#> 
#> pandoc 
#>   to: html
#>   output-file: report.html
#>   standalone: true
#>   section-divs: true
#>   html-math-method: mathjax
#>   wrap: none
#>   default-image-extension: png
#>   variables: {}
#>   
#> metadata
#>   document-css: false
#>   link-citations: true
#>   date-format: long
#>   lang: en
#>   engines:
#>     - path: /opt/quarto/share/extension-subtrees/julia-engine/_extensions/julia-engine/julia-engine.js
#>   title: MOSuite analysis report
#>   date: today
#>   
#> Output created: report.html
#> 
#> 
```
