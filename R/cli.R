# These functions were inspired by and adapted from renv:
#   https://github.com/rstudio/renv/blob/d0eb86349d35679eb6920ca59072bd7369fe620f/R/cli.R

#' Execute MOSuite from the CLI
#'
#' @export
#' @keywords internal
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
  # if --json in arguments, call cli_from_json()
  if (any(stringr::str_detect(clargs, "^--json"))) {
    args <- list(call(":::", as.symbol("MOSuite"), as.symbol("cli_from_json")),
      method = method
    )
  } else { # otherwise call the method directly
    args <- list(call("::", as.symbol("MOSuite"), as.symbol(method)))
  }

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

cli_usage <- function(con = stderr()) {
  usage <- "
Usage: mosuite [function] [--json=path/to/args.json]

[function] should be the name of a function exported from MOSuite.
[--json] should specify the path to a JSON file with arguments accepted by that function. The equals sign (=) is required to separate --json from the path.

Additionally, the JSON file can contain the following keys:
  - moo_input_rds: file path to an existing MultiOmicsDataset object in RDS format. This is required if `method` has `moo` as an argument.
  - moo_output_rds: file path to write the result to.

Use `mosuite [function] --help` for more information about the associated function.

Main functions:
  mosuite create_multiOmicDataSet_from_files
  mosuite filter_counts
  mosuite clean_raw_counts
  mosuite normalize_counts
  mosuite batch_correct_counts
  mosuite diff_counts
  mosuite filter_diff
"
  writeLines(usage, con = con)
}

cli_help <- function(method) {
  print(utils::help(method, package = "MOSuite"))
}

cli_unknown <- function(method, exports) {
  # report unknown command
  warning(glue::glue("MOSuite: {method} is not a known function."))

  # check for similar commands
  distance <- c(utils::adist(method, exports))
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

#' Call an MOSuite function with arguments specified in a json file
#'
#' @param method function in MOSuite to call
#' @param json path to a JSON file containing arguments for the function.
#'  Additionally, the JSON can contain the following keys:
#'    - `moo_input_rds` - filepath to an existing MultiOmicsDataset object in RDS format. This is required if the MOSuite function contains `moo` as an argument.
#'    - `moo_output_rds` - filepath to write the result to.
#' @param debug when TRUE, do not call the command, just return the expression.
#'
#' @returns invisible returns the function call
#' @export
#' @keywords internal
#'
cli_from_json <- function(method, json, debug = FALSE) {
  # begin building function call
  fcn_args <- list(call("::", as.symbol("MOSuite"), as.symbol(method)))
  # get function arguments from json
  json_args <- jsonlite::read_json(json)

  # if needed, get moo from moo_input_rds
  accepted_args <- formals(method, envir = getNamespace("MOSuite"))
  first_arg <- names(formals(method, envir = getNamespace("MOSuite")))[1]
  if (stringr::str_detect(first_arg, "^moo")) {
    assertthat::assert_that("moo_input_rds" %in% names(json_args),
      msg = glue::glue("moo_input_rds must be included in the JSON because `moo` is required for {method}()")
    )
    fcn_args[[first_arg]] <- readr::read_rds(json_args[["moo_input_rds"]])
  }
  # all other json keys should be arguments for the method
  # TODO convert lists to vectors
  fcn_args <- c(fcn_args, json_args[!stringr::str_detect(names(json_args), "moo_.*_rds")])

  # invoke method with parsed arguments from json
  expr <- as.call(fcn_args)
  if (isFALSE(debug)) {
    result <- eval(expr = expr, envir = globalenv())

    # save result to output_rds
    if ("moo_output_rds" %in% names(json_args)) {
      readr::write_rds(result, json_args[["moo_output_rds"]])
    }
  }

  invisible(expr)
}
