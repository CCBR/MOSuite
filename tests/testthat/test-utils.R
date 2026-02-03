test_that("abort_packages_not_installed works", {
  expect_no_error(abort_packages_not_installed("base"))
  expect_error(
    abort_packages_not_installed("not-a-package-name"),
    "The following package\\(s\\) are required but are not installed"
  )
})
test_that("check_packages_installed works", {
  expect_equal(
    check_packages_installed("base"),
    c(base = TRUE)
  )
  expect_equal(
    check_packages_installed("not-a-package-name"),
    c(`not-a-package-name` = FALSE)
  )
})
f <- function() stop("My error!")
expect_error(f())
expect_error(f(), "My error!")
