#!/usr/bin/env Rscript
# inherit the base image's repo (P3M binaries when available), else use CRAN
cran <- getOption('repos')[['CRAN']]
if (is.null(cran) || !startsWith(cran, 'http')) {
  options(repos = c(CRAN = 'https://cloud.r-project.org'))
}
pak::local_install_dev_deps('/opt/MOSuite', upgrade = FALSE)
dir.create('/data')
readr::write_csv(
  tibble::as_tibble(installed.packages()),
  '/data/r-packages_mosuite-dev.csv'
)
