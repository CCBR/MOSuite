#' @keywords internal
#' @examples
#'
#' get_function_meta(tools::Rd_db("MOSuite"), "batch_correct_counts")
#'
get_function_meta <- function(func_name, rd_db) {
  func_file <- paste0(func_name, ".Rd")
  func_db <- rd_db[[func_file]]

  title <- tools:::.Rd_get_metadata(func_db, "title")
  desc <- tools:::.Rd_get_metadata(func_db, "description")
  arg_desc <- dplyr::as_tibble(tools:::.Rd_get_argument_table(func_db), .name_repair = "unique_quiet")
  colnames(arg_desc) <- c("arg", "desc")
  arg_docs <- arg_desc %>%
    dplyr::pull("desc") %>%
    as.list()
  names(arg_docs) <- arg_desc %>% dplyr::pull("arg")
  # TODO Rd to markdown??

  arg_defaults <- lapply(formals(func_name), \(x) {
    if (inherits(x, "name")) {
      default <- NULL
    } else {
      default <- x
    }
    return(default)
  })
  if ("..." %in% names(arg_defaults)) {
    arg_defaults <- arg_defaults %>%
      within(., rm("...")) # remove `...` argument
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
#'   system.file("extdata", "galaxy", "template-templates", "create_multiOmicDataSet_from_files.json", package = "MOSuite"),
#'   tools::Rd_db("MOSuite")
#' )
#'
update_function_template <- function(template_filename, rd_db) {
  template <- jsonlite::read_json(template_filename)
  r_function <- template$r_function
  message(
    glue::glue(
      "Processing template {basename(template_filename)} for R function {r_function}"
    )
  )

  func_meta <- get_function_meta(r_function, rd_db)
  template$description <- func_meta$description
  args_in_template <- c()
  for (arg_type in c("columns", "inputDatasets", "parameters")) {
    for (i in seq_along(template[[arg_type]])) {
      arg_name <- template[[arg_type]][[i]]$key
      if (arg_name %in% names(func_meta$args)) {
        template[[arg_type]][[i]]$description <- func_meta$args[[arg_name]]$description
        template[[arg_type]][[i]]$defaultValue <- func_meta$args[[arg_name]]$defaultValue
        template[[arg_type]][[i]]$paramType <-
          args_in_template <- c(args_in_template, arg_name)
      } else {
        message(
          glue::glue(
            "{basename(template_filename)}: Argument {arg_name} from template not found in R function doc"
          )
        )
      }
    }
  }

  func_args_missing <- setdiff(names(func_meta$args), args_in_template)
  if (length(func_args_missing) > 0) {
    message(
      glue::glue(
        "{r_function}: Argument(s) from R function doc not found in template: {paste(func_args_missing, collapse = ', ')}"
      )
    )
  }
  return(template)
}

#' @keywords internal
check_classes <- function(updated_template) {
  for (p in updated_template$parameters) {
    for (el in p) {
      message(paste(p["key"], class(el)))
    }
  }
}

#' @keywords internal
write_package_json <- function(input_dir = file.path("inst", "extdata", "galaxy", "template-templates"),
                               output_dir = file.path("inst", "extdata", "galaxy", "code-templates")) {
  options(moo_print_plots = TRUE)
  options(moo_save_plots = TRUE)
  options(moo_plots_dir = "./figures")
  templates <- list.files(input_dir, pattern = ".*\\.json$", full.names = TRUE)
  rd_db <- tools::Rd_db("MOSuite")
  for (f in templates) {
    updated_template <- update_function_template(f, rd_db)
    jsonlite::write_json(
      updated_template,
      file.path(output_dir, basename(f)),
      auto_unbox = TRUE,
      pretty = TRUE,
      null = "null",
      na = "null"
    )
  }
}
