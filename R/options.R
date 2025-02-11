options::define_option(
  option = "print_plots",
  desc = "Whether to print plots during analysis",
  default = FALSE,
  option_name = "moo_print_plots",
  envvar_name = "MOO_PRINT_PLOTS",
  envvar_fn = options::envvar_is_true()
)

#' @eval options::as_roxygen_docs()
NULL

#' @eval options::as_params()
#' @name option_params
#'
NULL
