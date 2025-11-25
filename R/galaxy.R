#' @keywords internal
#' @examples
#'
#' get_function_meta(tools::Rd_db("MOSuite"), "batch_correct_counts")
#'
get_function_meta <- function(func_name, rd_db) {
  func_db <- rd_db[[paste0(func_name, ".Rd")]]

  title <- tools:::.Rd_get_metadata(func_db, "title")
  desc <- tools:::.Rd_get_metadata(func_db, "description")
  arg_desc <- dplyr::as_tibble(tools:::.Rd_get_argument_table(func_db), .name_repair = "unique_quiet")
  colnames(arg_desc) <- c("arg", "desc")
  arg_docs <- arg_desc %>%
    dplyr::pull("desc") %>%
    as.list()
  names(arg_docs) <- arg_desc %>% dplyr::pull("arg")

  arg_defaults <- lapply(formals(func_name), \(x) {
    if (inherits(x, "name")) {
      default <- NULL
    } else if (inherits(x, "call")) {
      default <- eval(x)
    } else {
      default <- x
    }
    return(default)
  })
  if ("..." %in% names(arg_defaults)) {
    arg_defaults <- arg_defaults %>%
      within(rm("...")) # remove `...` argument
  }
  args_meta <- names(arg_defaults) %>% lapply(\(arg) {
    return(list(defaultValue = arg_defaults[[arg]], description = arg_docs[[arg]]))
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
#' @examples
#'
#' update_function_template(
#'   system.file("extdata", "galaxy", "template-templates", "create_multiOmicDataSet_from_files.json",
#'     package = "MOSuite"
#'   ),
#'   tools::Rd_db("MOSuite")
#' )
#'
update_function_template <- function(template,
                                     func_meta,
                                     func_defaults,
                                     keep_deprecated_args = TRUE) {
  abort_packages_not_installed("Rd2md")
  new_template <- list(
    r_function = template$r_function,
    title = template$title |> Rd2md::rd_str_to_md(),
    description = func_meta$description |> Rd2md::rd_str_to_md(),
    columns = list(),
    inputDatasets = list(),
    parameters = list()
  )
  args_in_template <- c()
  template_args_missing <- c()
  for (arg_type in c("columns", "inputDatasets", "parameters")) {
    for (i in seq_along(template[[arg_type]])) {
      arg_name <- template[[arg_type]][[i]]$key
      if (arg_name %in% names(func_meta$args)) {
        arg_meta <- template[[arg_type]][[i]]
        arg_meta$description <- func_meta$args[[arg_name]]$description |> Rd2md::rd_str_to_md()
        arg_meta$defaultValue <- func_defaults[[arg_name]]
        args_in_template <- c(args_in_template, arg_name)
        new_template[[arg_type]][[length(new_template[[arg_type]]) + 1]] <- arg_meta
      } else {
        template_args_missing <- c(template_args_missing, arg_name)
        if (isTRUE(keep_deprecated_args)) {
          arg_meta <- template[[arg_type]][[i]]
          new_template[[arg_type]][[length(new_template[[arg_type]]) + 1]] <- arg_meta
        }
      }
    }
  }

  if (length(template_args_missing) > 0) {
    glue::glue(
      "{basename(template_filename)}: ",
      "Argument(s) from template not found in R function doc: ",
      "{paste(template_args_missing, collapse = ', ')}"
    )
  }

  func_args_missing <- setdiff(names(func_meta$args), args_in_template)
  if (length(func_args_missing) > 0) {
    message(
      glue::glue(
        "{r_function}: ",
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
write_json <- function(x,
                       filepath,
                       auto_unbox = TRUE,
                       pretty = TRUE,
                       null = "null",
                       na = "null",
                       ...) {
  jsonlite::write_json(
    x,
    filepath,
    auto_unbox = auto_unbox,
    pretty = pretty,
    null = null,
    na = na,
    ...
  )
}

#' @keywords internal
get_function_defaults <- function(func_meta) {
  func_defaults <- list()
  names(func_meta$args) %>% purrr::map(\(arg_name) {
    if (stringr::str_starts(arg_name, "moo")) {
      func_defaults[["moo_input_rds"]] <- "moo.rds"
      func_defaults[["moo_output_rds"]] <- "moo.rds"
    } else {
      func_defaults[[arg_name]] <- func_meta$args[[arg_name]]$defaultValue
    }
  })

  return(func_defaults)
}

#' @keywords internal
write_package_json_blueprints <- function(input_dir = file.path("inst", "extdata", "galaxy", "template-templates"),
                                          blueprints_output_dir = file.path("inst", "extdata", "galaxy", "code-templates"),
                                          defaults_output_dir = file.path("inst", "extdata", "json_args", "defaults")) {
  options(moo_print_plots = TRUE)
  options(moo_save_plots = TRUE)
  options(moo_plots_dir = "./figures")
  templates <- list.files(input_dir, pattern = ".*\\.json$", full.names = TRUE)
  rd_db <- tools::Rd_db("MOSuite")
  for (f in templates) {
    base_filename <- basename(f)
    message(glue::glue("* Processing {base_filename}"))

    template <- jsonlite::read_json(f)
    r_function <- template$r_function
    func_meta <- get_function_meta(r_function, rd_db)

    # write default arguments
    func_defaults <- get_function_defaults(func_meta)
    write_json(
      func_defaults,
      file.path(defaults_output_dir, glue::glue("{r_function}.json"))
    )

    # write galaxy blueprint template
    updated_template <- update_function_template(template, func_meta, func_defaults)
    write_json(
      updated_template,
      file.path(blueprints_output_dir, glue::glue("{r_function}.json"))
    )
  }
  return()
}
