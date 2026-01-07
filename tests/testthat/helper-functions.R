equal_dfs <- function(x, y) {
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
      tmp <- capture.output(
        test_that(desc = desc, ...)
      )
    }
  )
  return(invisible())
}
