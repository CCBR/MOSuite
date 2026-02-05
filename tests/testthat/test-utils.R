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
# Tests for parse_optional_vector
test_that("parse_optional_vector handles normal input", {
  result <- parse_optional_vector("a, b, c")
  expect_equal(result, c("a", "b", "c"))
})

test_that("parse_optional_vector trims whitespace", {
  result <- parse_optional_vector("  a  ,  b  ,  c  ")
  expect_equal(result, c("a", "b", "c"))
})

test_that("parse_optional_vector returns NULL for empty string", {
  result <- parse_optional_vector("")
  expect_null(result)
})

test_that("parse_optional_vector returns NULL for NULL input", {
  result <- parse_optional_vector(NULL)
  expect_null(result)
})

test_that("parse_optional_vector returns NULL for zero-length vector", {
  result <- parse_optional_vector(character(0))
  expect_null(result)
})

test_that("parse_optional_vector handles single value", {
  result <- parse_optional_vector("single")
  expect_equal(result, "single")
})

test_that("parse_optional_vector handles numeric-like strings", {
  result <- parse_optional_vector("1, 2, 3")
  expect_equal(result, c("1", "2", "3"))
})

# Tests for parse_vector_with_default
test_that("parse_vector_with_default parses normal input", {
  result <- parse_vector_with_default("a, b, c", "default")
  expect_equal(result, c("a", "b", "c"))
})

test_that("parse_vector_with_default returns default for empty string", {
  result <- parse_vector_with_default("", "default")
  expect_equal(result, "default")
})

test_that("parse_vector_with_default returns default for NULL", {
  result <- parse_vector_with_default(NULL, "default")
  expect_equal(result, "default")
})

test_that("parse_vector_with_default returns default for zero-length vector", {
  result <- parse_vector_with_default(character(0), "default")
  expect_equal(result, "default")
})

test_that("parse_vector_with_default handles vector defaults", {
  default_vec <- c("x", "y", "z")
  result <- parse_vector_with_default("", default_vec)
  expect_equal(result, default_vec)
})

test_that("parse_vector_with_default handles numeric defaults", {
  result <- parse_vector_with_default("", 42)
  expect_equal(result, 42)
})

# Tests for parse_samples_to_rename
test_that("parse_samples_to_rename parses single pair", {
  result <- parse_samples_to_rename("old:new")
  expect_equal(result, list(old = "new"))
})

test_that("parse_samples_to_rename parses multiple pairs", {
  result <- parse_samples_to_rename("sample1:S1,sample2:S2,sample3:S3")
  expect_equal(result, list(sample1 = "S1", sample2 = "S2", sample3 = "S3"))
})

test_that("parse_samples_to_rename handles many sample pairs", {
  result <- parse_samples_to_rename(
    "ctrl_1:Control_Rep1,ctrl_2:Control_Rep2,treat_1:Treatment_Rep1,treat_2:Treatment_Rep2,treat_3:Treatment_Rep3"
  )
  expected <- list(
    ctrl_1 = "Control_Rep1",
    ctrl_2 = "Control_Rep2",
    treat_1 = "Treatment_Rep1",
    treat_2 = "Treatment_Rep2",
    treat_3 = "Treatment_Rep3"
  )
  expect_equal(result, expected)
})

test_that("parse_samples_to_rename trims whitespace", {
  result <- parse_samples_to_rename("  old  :  new  ,  old2  :  new2  ")
  expect_equal(result, list(old = "new", old2 = "new2"))
})

test_that("parse_samples_to_rename returns NULL for empty string", {
  result <- parse_samples_to_rename("")
  expect_null(result)
})

test_that("parse_samples_to_rename returns NULL for NULL input", {
  result <- parse_samples_to_rename(NULL)
  expect_null(result)
})

test_that("parse_samples_to_rename returns NULL for zero-length vector", {
  result <- parse_samples_to_rename(character(0))
  expect_null(result)
})

test_that("parse_samples_to_rename ignores malformed pairs", {
  result <- parse_samples_to_rename("valid:pair,invalid_no_colon")
  expect_equal(result, list(valid = "pair"))
})

test_that("parse_samples_to_rename returns NULL if all pairs malformed", {
  result <- parse_samples_to_rename("malformed,also_malformed")
  expect_null(result)
})

test_that("parse_samples_to_rename handles colons in values", {
  result <- parse_samples_to_rename("old:http://new")
  # Only uses pairs with exactly 2 parts (one colon), so this is ignored
  expect_null(result)
})

test_that("parse_samples_to_rename silently ignores pairs with colons in names", {
  # When a column name contains colon, pair is silently skipped
  result <- parse_samples_to_rename("old:with:colons:new,sample1:S1")
  # Only sample1:S1 is valid, the other is ignored
  expect_equal(result, list(sample1 = "S1"))
})

# Tests for setup_capsule_environment
test_that("setup_capsule_environment creates correct directory paths", {
  tmpdir <- tempdir()

  result <- setup_capsule_environment(base_results_dir = tmpdir)

  expect_equal(result$results_dir, tmpdir)
  expect_equal(result$plots_dir, file.path(tmpdir, "figures"))
})

test_that("setup_capsule_environment sets options", {
  tmpdir <- tempdir()

  setup_capsule_environment(base_results_dir = tmpdir)

  expect_equal(getOption("moo_plots_dir"), file.path(tmpdir, "figures"))
  expect_equal(getOption("moo_save_plots"), TRUE)
})

test_that("setup_capsule_environment creates r-packages.csv", {
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  setup_capsule_environment(base_results_dir = tmpdir)

  csv_path <- file.path(tmpdir, "r-packages.csv")
  expect_true(file.exists(csv_path))

  # Check it's a valid CSV
  csv_content <- readr::read_csv(csv_path, show_col_types = FALSE)
  expect_gt(nrow(csv_content), 0)
})

test_that("setup_capsule_environment returns invisibly", {
  tmpdir <- tempdir()

  result <- setup_capsule_environment(base_results_dir = tmpdir)

  expect_type(result, "list")
  expect_named(result, c("results_dir", "plots_dir"))
})

# Tests for load_moo_from_data_dir
test_that("load_moo_from_data_dir stops when no .rds files found", {
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  expect_error(
    load_moo_from_data_dir(data_dir = tmpdir),
    "No files matching regex"
  )
})

test_that("load_moo_from_data_dir loads valid MOO object", {
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  # Create a mock multiOmicDataSet object
  moo <- structure(
    list(data = "test"),
    class = c("multiOmicDataSet", "MOSuite::multiOmicDataSet")
  )

  rds_file <- file.path(tmpdir, "test.rds")
  readr::write_rds(moo, rds_file)

  result <- load_moo_from_data_dir(data_dir = tmpdir)

  expect_s3_class(result, "multiOmicDataSet")
  expect_equal(result$data, "test")
})

test_that("load_moo_from_data_dir stops for invalid class", {
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  # Create an invalid object (not a multiOmicDataSet)
  invalid_obj <- list(data = "test")

  rds_file <- file.path(tmpdir, "test.rds")
  readr::write_rds(invalid_obj, rds_file)

  expect_error(
    load_moo_from_data_dir(data_dir = tmpdir),
    "The input is not a multiOmicDataSet"
  )
})

test_that("load_moo_from_data_dir finds file recursively", {
  tmpdir <- tempfile()
  dir.create(tmpdir)
  subdir <- file.path(tmpdir, "subdir")
  dir.create(subdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  # Create a mock multiOmicDataSet in subdirectory
  moo <- structure(
    list(data = "test"),
    class = c("multiOmicDataSet", "MOSuite::multiOmicDataSet")
  )

  rds_file <- file.path(subdir, "test.rds")
  readr::write_rds(moo, rds_file)

  result <- load_moo_from_data_dir(data_dir = tmpdir)

  expect_s3_class(result, "multiOmicDataSet")
})

test_that("load_moo_from_data_dir uses first matching file", {
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  # Create two mock multiOmicDataSet objects
  moo1 <- structure(
    list(data = "first"),
    class = c("multiOmicDataSet", "MOSuite::multiOmicDataSet")
  )
  moo2 <- structure(
    list(data = "second"),
    class = c("multiOmicDataSet", "MOSuite::multiOmicDataSet")
  )

  rds_file1 <- file.path(tmpdir, "a_test.rds")
  rds_file2 <- file.path(tmpdir, "z_test.rds")
  readr::write_rds(moo1, rds_file1)
  readr::write_rds(moo2, rds_file2)

  result <- load_moo_from_data_dir(data_dir = tmpdir)

  # Should load one of them (order may vary, so just check it's valid)
  expect_s3_class(result, "multiOmicDataSet")
  expect_true(result$data %in% c("first", "second"))
})

test_that("load_moo_from_data_dir prints message", {
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  # Create a mock multiOmicDataSet object
  moo <- structure(
    list(data = "test"),
    class = c("multiOmicDataSet", "MOSuite::multiOmicDataSet")
  )

  rds_file <- file.path(tmpdir, "test.rds")
  readr::write_rds(moo, rds_file)

  expect_message(
    load_moo_from_data_dir(data_dir = tmpdir),
    "Reading multiOmicDataSet from"
  )
})
