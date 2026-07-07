# source: https://rconsortium.github.io/S7/articles/packages.html#method-registration
.onLoad <- function(...) {
  S7::methods_register()

  # necessary because moo_save_plots defaults to TRUE, causing figures to be
  # written during R CMD check
  is_r_cmd_check <- nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")) ||
    identical(tolower(Sys.getenv("R_CMD_CHECK")), "true") ||
    identical(tolower(Sys.getenv("R CMD check")), "true")
  if (
    is_r_cmd_check &&
      is.null(getOption("moo_save_plots")) &&
      !nzchar(Sys.getenv("MOO_SAVE_PLOTS"))
  ) {
    options(moo_save_plots = FALSE)
  }
}

# enable usage of <S7_object>@name in package code
# source: https://rconsortium.github.io/S7/articles/packages.html#backward-compatibility
#' @rawNamespace if (getRversion() < "4.3.0") importFrom("S7", "@")
NULL


# Suppress R CMD check note 'All declared Imports should be used'.
# These packages are used within S7 methods.
#' @importFrom DESeq2 DESeq
NULL
