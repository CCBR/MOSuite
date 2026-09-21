#!/usr/bin/env Rscript
options(repos = c(CRAN = 'https://cloud.r-project.org'))
pak::local_install_dev_deps('/opt/MOSuite', upgrade = FALSE)
dir.create('/data')
readr::write_csv(
  tibble::as_tibble(installed.packages()),
  '/data/r-packages_mosuite-dev.csv'
)
