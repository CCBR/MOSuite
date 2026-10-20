#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (
  length(arguments) != 2L ||
    !arguments[[1L]] %in% c("mosuite", "mosuite-minimal")
) {
  stop("Usage: check-image-packages.R <mosuite|mosuite-minimal> <DESCRIPTION>")
}

image_kind <- arguments[[1L]]
description <- read.dcf(
  arguments[[2L]],
  fields = c("Depends", "Imports", "LinkingTo")
)
dependency_fields <- description[1L, ]
dependency_fields <- dependency_fields[!is.na(dependency_fields)]
dependency_specs <- trimws(unlist(strsplit(
  paste(dependency_fields, collapse = ","),
  ",",
  fixed = TRUE
)))
dependency_names <- unique(trimws(sub(
  "[[:space:]]*[(].*$",
  "",
  dependency_specs
)))
dependency_names <- setdiff(dependency_names, c("R", "MOSuite"))
missing_dependencies <- dependency_names[
  !vapply(
    dependency_names,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]
if (length(missing_dependencies) > 0L) {
  stop(
    "Missing DESCRIPTION dependencies: ",
    paste(missing_dependencies, collapse = ", ")
  )
}

has_mosuite <- requireNamespace("MOSuite", quietly = TRUE)
should_have_mosuite <- identical(image_kind, "mosuite")
if (!identical(has_mosuite, should_have_mosuite)) {
  stop("Unexpected MOSuite installation for image: ", image_kind)
}
