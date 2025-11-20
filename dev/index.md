# MOSuite

R package for differential multi-omics analysis

Multi-Omics Suite provides a suite of functions to clean, filter,
batch-correct, normalize, visualize, and perform differential expression
analysis. While the package is designed for differential
[RNA-seq](https://github.com/CCBR/RENEE) analysis and multi-omics
datasets, it can be used for any data represented in a counts table. See
the website for more information, documentation, and examples:
<https://ccbr.github.io/MOSuite/>

## Installation

You can install the development version of MOSuite from
[GitHub](https://github.com/CCBR/MOSuite) with:

``` r
# install.packages("remotes")
remotes::install_github("CCBR/MOSuite", dependencies = TRUE)
```

Or install a specific version:

``` r
remotes::install_github("CCBR/MOSuite", dependencies = TRUE, ref = "v0.1.0")
```

There is also a Docker container available at
<https://hub.docker.com/r/nciccbr/mosuite>

``` sh
singularity exec docker://nciccbr/mosuite:v0.2.0 R -s -e \
  'cat("MOSuite version:", installed.packages()["MOSuite",][["Version"]])'
```

## Usage

Please see the [introductory
vignette](https://ccbr.github.io/MOSuite/articles/intro.html) for a
quick start tutorial. Take a look at the [reference
documentation](https://ccbr.github.io/MOSuite/reference/index.html) for
detailed information on each function in the package.

## Help & Contributing

Come across a **bug**? Open an
[issue](https://github.com/CCBR/MOSuite/issues) and include a minimal
reproducible example.

Have a **question**? Ask it in
[discussions](https://github.com/CCBR/MOSuite/discussions).

Want to **contribute** to this project? Check out the [contributing
guidelines](https://ccbr.github.io/MOSuite/dev/CONTRIBUTING.md).

## Development Roadmap

![](./reference/figures/development-plan.png)

- [dev
  spreadsheet](https://nih-my.sharepoint.com/:x:/g/personal/homanpj_nih_gov/ETvHXgnwxExEpcP57Jj9_EwBHBvZBqNuZ_c3eu51w-SlnA?e=PcXKU8)
- [project board](https://github.com/orgs/CCBR/projects/32)
