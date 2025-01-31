equal_dfs <- function(x, y) {
  all(
    class(x) == class(y),
    names(x) == names(y),
    rownames(x) == rownames(y),
    all.equal(x, y),
    all.equal(lapply(x, class), lapply(y, class))
  )
}
