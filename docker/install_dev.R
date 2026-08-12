#!/usr/bin/env Rscript
options(repos = c(CRAN = 'https://cloud.r-project.org'))
install.packages('pak')
pak::local_install_dev_deps('/opt/MOSuite', upgrade = TRUE)
dir.create('/data')
readr::write_csv(
  tibble::as_tibble(installed.packages()),
  '/data/r-packages_mosuite-minimal.csv'
)
