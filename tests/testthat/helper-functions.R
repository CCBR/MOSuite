equal_dfs <- function(x, y) {
  all(
    class(x) == class(y),
    names(x) == names(y),
    rownames(x) == rownames(y),
    all.equal(x, y),
    all.equal(lapply(x, class), lapply(y, class))
  )
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

  list(object = x, path = paste0("compare_proxy(", path, ")"))
}
