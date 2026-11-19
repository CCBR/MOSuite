equal_dfs <- function(x, y) {
  # Ignore readr parse metadata so comparisons focus on data content.
  attr(x, "spec") <- NULL
  attr(x, "problems") <- NULL
  attr(y, "spec") <- NULL
  attr(y, "problems") <- NULL

  return(all(
    class(x) == class(y),
    names(x) == names(y),
    rownames(x) == rownames(y),
    all.equal(x, y),
    all.equal(lapply(x, class), lapply(y, class))
  ))
}

# source https://stackoverflow.com/a/75232781/5787827
compare_proxy.plotly <- function(x, path = "x") {
  names(x$x$visdat) <- "proxy"
  e <- environment(x$x$visdat$proxy)

  # Maybe we should follow the recursion, but not now.
  e$p <- NULL

  e$id <- "proxy"

  x$x$cur_data <- "proxy"
  names(x$x$attrs) <- "proxy"

  return(list(object = x, path = paste0("compare_proxy(", path, ")")))
}

run_function_cli <- function(func_name) {
  json_path <- paste0(
    func_name,
    ".json"
  )

  return(cli_exec(c(
    func_name,
    paste0('--json="', json_path, '"')
  )))
}

# source: https://github.com/r-lib/testthat/issues/664#issuecomment-340809997
create_empty_dir <- function(x) {
  unlink(x, recursive = TRUE, force = TRUE)
  return(dir.create(x))
}

# source: https://github.com/r-lib/testthat/issues/664#issuecomment-340809997
test_with_dir <- function(desc, ...) {
  new <- tempfile()
  create_empty_dir(new)
  withr::with_dir(
    # or local_dir()
    new = new,
    code = {
      capture.output(
        testthat::test_that(desc = desc, ...) # nolint: object_usage_linter
      )
    }
  )
  return(invisible())
}

pca_point_coordinates <- function(plot) {
  built_plot <- ggplot2::ggplot_build(plot)
  return(built_plot$data[[1]][, c("x", "y")])
}

capture_saved_pca_plot <- function() {
  pca_plot <- NULL
  print_or_save_plot <- function(plot, filename, ...) {
    if (basename(filename) == "pca.png") {
      pca_plot <<- plot
    }
    return(invisible(NULL))
  }
  return(list(
    get = function() pca_plot,
    print_or_save_plot = print_or_save_plot
  ))
}

expect_pca_coordinates_equal <- function(actual_plot, expected_plot) {
  expect_equal(
    pca_point_coordinates(actual_plot),
    pca_point_coordinates(expected_plot),
    tolerance = 1e-8
  )
}

histogram_layer_data <- function(plot) {
  built_plot <- ggplot2::ggplot_build(plot)
  return(built_plot$data)
}

capture_saved_histogram_plot <- function() {
  histogram_plot <- NULL
  print_or_save_plot <- function(plot, filename, ...) {
    if (basename(filename) == "histogram.png") {
      histogram_plot <<- plot
    }
    return(invisible(NULL))
  }
  return(list(
    get = function() histogram_plot,
    print_or_save_plot = print_or_save_plot
  ))
}

expect_histogram_layers_equal <- function(actual_plot, expected_plot) {
  expect_equal(
    histogram_layer_data(actual_plot),
    histogram_layer_data(expected_plot),
    tolerance = 1e-8
  )
}

create_tiny_moo <- function(tag = "test") {
  sample_metadata <- data.frame(sample_id = c("s1", "s2"), group = c("A", "B"))
  counts_dat <- data.frame(
    feature_id = c("g1", "g2"),
    s1 = c(1, 2),
    s2 = c(3, 4)
  )
  anno_dat <- data.frame(feature_id = c("g1", "g2"), symbol = c("G1", "G2"))

  moo <- multiOmicDataSet(
    sample_metadata = sample_metadata,
    anno_dat = anno_dat,
    counts_lst = list(raw = counts_dat)
  )
  moo@analyses$tag <- tag

  return(moo)
}
