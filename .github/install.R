#!/usr/bin/env Rscript --vanilla
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("remotes")
remotes::install_version("ggplot2", version = "3.5.2")
