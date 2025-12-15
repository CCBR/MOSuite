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
  json_path <- testthat::test_path("data", paste0(
    func_name, ".json"
  ))
  return(cli_exec(c(
    func_name, paste0('--json="', json_path, '"')
  )))
}
