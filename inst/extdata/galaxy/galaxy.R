library(dplyr)

#' @keywords internal
#' @examples
#'
#' get_function_meta("batch_correct_counts", tools::Rd_db("MOSuite"))
#'
get_function_meta <- function(func_name, rd_db) {
  func_db <- rd_db[[paste0(func_name, ".Rd")]]

  title <- tools:::.Rd_get_metadata(func_db, "title") |> trimws()
  desc <- paste(
    tools:::.Rd_get_metadata(func_db, "description"),
    tools:::.Rd_get_metadata(func_db, "details"),
    sep = "\n\n"
  ) |>
    trimws()
  arg_desc <- dplyr::as_tibble(
    tools:::.Rd_get_argument_table(func_db),
    .name_repair = "unique_quiet"
  )
  colnames(arg_desc) <- c("arg", "desc")
  arg_docs <- arg_desc |>
    dplyr::pull("desc") |>
    trimws() |>
    as.list()
  names(arg_docs) <- arg_desc |> dplyr::pull("arg")
  options(
    moo_print_plots = TRUE,
    moo_save_plots = TRUE,
    moo_plots_dir = "./figures",
    print_plots = TRUE, # need if this function is defined outside the package's R source directory
    save_plots = TRUE,
    plots_dir = "./figures"
  )
  arg_defaults <- lapply(
    formals(func_name, envir = getNamespace("MOSuite")),
    \(x) {
      if (inherits(x, "name")) {
        default <- NULL
      } else if (inherits(x, "call")) {
        default <- eval(x, envir = getNamespace("MOSuite"))
      } else {
        default <- x
      }
      return(default)
    }
  )
  if ("..." %in% names(arg_defaults)) {
    arg_defaults <- arg_defaults |>
      within(rm("...")) # remove `...` argument
  }
  args_meta <- names(arg_defaults) |>
    lapply(\(arg) {
      return(list(
        defaultValue = arg_defaults[[arg]],
        description = arg_docs[[arg]]
      ))
    })
  names(args_meta) <- names(arg_defaults)

  return(list(
    r_function = func_name,
    title = title,
    description = desc,
    args = args_meta
  ))
}

#' @keywords internal
get_function_args <- function(func_meta) {
  func_names <- Filter(
    \(x) !stringr::str_starts(x, "moo"),
    names(func_meta$args)
  )
  func_args <- lapply(func_names, \(x) func_meta$args[[x]][["defaultValue"]])

  if (stringr::str_starts(names(func_meta$args)[1], "moo")) {
    func_names <- c("moo_input_rds", "moo_output_rds", func_names)
    func_args <- c("moo.rds", "moo.rds", func_args)
  }
  names(func_args) <- func_names

  return(func_args)
}

#' @keywords internal
#' @examples
#'
#' update_function_template(
#'   system.file("extdata", "galaxy", "template-templates", "create_multiOmicDataSet_from_files.json",
#'     package = "MOSuite"
#'   ),
#'   tools::Rd_db("MOSuite")
#' )
#'
update_function_template <- function(
  template,
  func_meta,
  keep_deprecated_args = TRUE
) {
  if (!rlang::is_installed("Rd2md")) {
    stop("Required pacakge {Rd2md} is not installed")
  }
  new_template <- list(
    r_function = template$r_function,
    title = template$title |> Rd2md::rd_str_to_md(),
    description = func_meta$description |> Rd2md::rd_str_to_md(),
    columns = list(),
    inputDatasets = list(),
    parameters = list(),
    outputs = template$outputs
  )
  args_in_template <- c()
  template_args_missing <- c()
  for (arg_type in c("columns", "inputDatasets", "parameters")) {
    for (i in seq_along(template[[arg_type]])) {
      arg_name <- template[[arg_type]][[i]]$key
      if (arg_name %in% names(func_meta$args)) {
        arg_meta <- template[[arg_type]][[i]]
        arg_meta$description <- func_meta$args[[arg_name]]$description |>
          Rd2md::rd_str_to_md()
        arg_meta$defaultValue <- func_meta$args[[arg_name]]$defaultValue
        args_in_template <- c(args_in_template, arg_name)
        new_template[[arg_type]][[
          length(new_template[[arg_type]]) + 1
        ]] <- arg_meta
      } else {
        template_args_missing <- c(template_args_missing, arg_name)
        if (isTRUE(keep_deprecated_args)) {
          arg_meta <- template[[arg_type]][[i]]
          new_template[[arg_type]][[
            length(new_template[[arg_type]]) + 1
          ]] <- arg_meta
        }
      }
    }
  }
  if (length(template_args_missing) > 0) {
    message(glue::glue(
      "{template$r_function}: ",
      "Argument(s) from template not found in R function doc: ",
      "{paste(template_args_missing, collapse = ', ')}"
    ))
  }

  func_args_missing <- setdiff(names(func_meta$args), args_in_template)
  if (length(func_args_missing) > 0) {
    message(
      glue::glue(
        "{template$r_function}: ",
        "Argument(s) from R function doc not found in template: ",
        "{paste(func_args_missing, collapse = ', ')}"
      )
    )
  }
  return(new_template)
}

#' @keywords internal
check_classes <- function(updated_template) {
  for (p in updated_template$parameters) {
    for (el in p) {
      message(paste(p["key"], class(el)))
    }
  }
  return()
}

#' `jsonlite::write_json()` with preferred defaults
#'
#' @keywords internal
write_json <- function(
  x,
  filepath,
  auto_unbox = TRUE,
  pretty = TRUE,
  null = "null",
  na = "null",
  ...
) {
  return(invisible(jsonlite::write_json(
    x,
    filepath,
    auto_unbox = auto_unbox,
    pretty = pretty,
    null = null,
    na = na,
    ...
  )))
}

#' @keywords internal
write_package_json_blueprints <-
  function(
    input_dir = file.path("inst", "extdata", "galaxy", "1_mosuite-templates"),
    blueprints_output_dir = file.path(
      "inst",
      "extdata",
      "galaxy",
      "2_blueprints"
    ),
    defaults_output_dir = file.path("inst", "extdata", "json_args", "defaults")
  ) {
    options(
      moo_print_plots = TRUE,
      moo_save_plots = TRUE,
      moo_plots_dir = "./figures",
      print_plots = TRUE, # need if this function is defined outside the package's R source directory
      save_plots = TRUE,
      plots_dir = "./figures"
    )
    templates <- list.files(
      input_dir,
      pattern = ".*\\.json$",
      full.names = TRUE
    )
    rd_db <- tools::Rd_db("MOSuite")
    for (f in templates) {
      base_filename <- basename(f)
      message(glue::glue("* Processing {base_filename}"))

      template <- jsonlite::read_json(f)
      r_function <- template$r_function
      func_meta <- get_function_meta(r_function, rd_db)

      # write default arguments
      func_args <- get_function_args(func_meta)
      write_json(
        func_args,
        file.path(defaults_output_dir, glue::glue("{r_function}.json"))
      )

      # write galaxy blueprint template
      updated_template <- update_function_template(
        template,
        func_meta
      )
      write_json(
        updated_template,
        file.path(blueprints_output_dir, glue::glue("{r_function}.json"))
      )
    }
    return(invisible())
  }
