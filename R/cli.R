# These functions were inspired by and adapted from renv:
#   https://github.com/rstudio/MOSuite/blob/d0eb86349d35679eb6920ca59072bd7369fe620f/R/cli.R

#' Execute MOSuite from the CLI
#'
cli_exec <- function(clargs = commandArgs(trailingOnly = TRUE)) {
  invisible(cli_exec_impl(clargs))
}
cli_exec_impl <- function(clargs) {
  # check for tool called without arguments, or called with '--help'
  usage <-
    length(clargs) == 0 ||
      clargs[1L] %in% c("help", "--help")

  if (usage) {
    return(cli_usage())
  }

  # extract method
  method <- clargs[1L]

  # check request for help on requested method
  help <-
    clargs[2L] %in% c("help", "--help")

  if (help) {
    return(cli_help(method))
  }

  # check for known function in MOSuite
  exports <- getNamespaceExports("MOSuite")
  if (!method %in% exports) {
    return(cli_unknown(method, exports))
  }

  # begin building call
  args <- list(call("::", as.symbol("MOSuite"), as.symbol(method)))

  for (clarg in clargs[-1L]) {
    # convert '--no-<flag>' into a FALSE parameter
    if (grepl("^--no-", clarg)) {
      key <- substring(clarg, 6L)
      args[[key]] <- FALSE
    }

    # convert '--param=value' flags
    else if (grepl("^--[^=]+=", clarg)) {
      index <- regexpr("=", clarg, fixed = TRUE)
      key <- substring(clarg, 3L, index - 1L)
      val <- substring(clarg, index + 1L)
      args[[key]] <- cli_parse(val)
    }

    # convert '--flag' into a TRUE parameter
    else if (grepl("^--", clarg)) {
      key <- substring(clarg, 3L)
      args[[key]] <- TRUE
    }

    # convert 'param=value' flags
    else if (grepl("=", clarg, fixed = TRUE)) {
      index <- regexpr("=", clarg, fixed = TRUE)
      key <- substring(clarg, 1L, index - 1L)
      val <- substring(clarg, index + 1L)
      args[[key]] <- cli_parse(val)
    }

    # take other parameters as-is
    else {
      args[[length(args) + 1L]] <- cli_parse(clarg)
    }
  }

  # invoke method with parsed arguments
  expr <- as.call(args)
  eval(expr = expr, envir = globalenv())
}

cli_usage <- function() {
  usage <- "
Usage: mosuite [method] [args...]

[method] should be the name of a function exported from MOSuite.
[args...] should be arguments accepted by that function.

Use mosuite [method] --help for more information about the associated function.

Examples:

  # basic commands
  mosuite create_multiOmicDataSet_from_files   # Create MO object from files
"

  writeLines(usage, con = stderr())
}
cli_help <- function(method) {
  print(help(method, package = "MOSuite"))
}

cli_unknown <- function(method, exports) {
  # report unknown command
  warning("MOSuite: '%s' is not a known command.", method)

  # check for similar commands
  distance <- c(adist(method, exports))
  names(distance) <- exports
  n <- min(distance)
  if (n > 2) {
    return(1L)
  }

  candidates <- names(distance)[distance == n]
  fmt <- "did you mean %s?"
  warning(fmt, paste(shQuote(candidates), collapse = " or "))
  return(1L)
}

cli_parse <- function(text) {
  # handle logical-like values up-front
  if (text %in% c("true", "True", "TRUE")) {
    return(TRUE)
  } else if (text %in% c("false", "False", "FALSE")) {
    return(FALSE)
  }

  # parse the expression
  value <- parse(text = text)[[1L]]
  if (is.language(value)) text else value
}
