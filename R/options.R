options::set_option_name_fn(function(package, name) {
  tolower(paste0("moo_", name))
})

options::set_envvar_name_fn(function(package, name) {
  gsub("[^A-Z0-9]", "_", toupper(paste0("moo_", name)))
})

options::define_option(
  option = "print_plots",
  default = FALSE,
  desc = "Whether to print plots during analysis",
  option_name = "moo_print_plots",
  envvar_name = "MOO_PRINT_PLOTS"
)


#' @eval options::as_roxygen_docs()
NULL

#' @eval options::as_params()
#' @name option_params
#'
NULL
