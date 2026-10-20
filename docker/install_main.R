#!/usr/bin/env Rscript

cran <- getOption("repos")[["CRAN"]]
if (!is.null(cran) && grepl("__linux__", cran, fixed = TRUE)) {
  options(repos = c(CRAN = sub("/[^/]+$", "/latest", cran)))
} else {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

pak::pkg_install(c("devtools", "readr", "tibble"))
pak::pkg_install("local::/opt/MOSuite", dependencies = TRUE, upgrade = FALSE)
devtools::test("/opt/MOSuite", stop_on_failure = TRUE)

dir.create("/data", showWarnings = FALSE)
readr::write_csv(
  tibble::as_tibble(installed.packages()),
  "/data/r-packages_mosuite.csv"
)
