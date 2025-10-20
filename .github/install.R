#!/usr/bin/env Rscript --vanilla
install.packages("remotes")
remotes::install_version("ggplot2", version = "3.5.2")
remotes::install_local(dependencies = TRUE, upgrade = "never")
