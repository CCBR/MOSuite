#!/usr/bin/env Rscript
# rocker pins a dated P3M snapshot, which is too old for some dependencies.
# Use the latest P3M snapshot to keep binary packages, or CRAN if the base
# image has no P3M repo configured (e.g. arm64, which P3M does not serve).
cran <- getOption('repos')[['CRAN']]
if (!is.null(cran) && grepl('__linux__', cran, fixed = TRUE)) {
  options(repos = c(CRAN = sub('/[^/]+$', '/latest', cran)))
} else {
  options(repos = c(CRAN = 'https://cloud.r-project.org'))
}
pak::local_install_dev_deps('/opt/MOSuite', upgrade = FALSE)
dir.create('/data')
readr::write_csv(
  tibble::as_tibble(installed.packages()),
  '/data/r-packages_mosuite-dev.csv'
)
