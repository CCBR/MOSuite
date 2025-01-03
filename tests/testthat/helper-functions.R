equal_dfs <- function(x, y) {
  all(
    class(x) == class(y),
    names(x) == names(y),
    x == y,
    all.equal(lapply(x, class), lapply(y, class))
  )
}
