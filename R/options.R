options::set_option_name_fn(function(package, name) {
  return(tolower(paste0("moo_", name)))
})

options::set_envvar_name_fn(function(package, name) {
  return(gsub("[^A-Z0-9]", "_", toupper(paste0("moo_", name))))
})

options::define_option(
  option = "print_plots",
  default = FALSE,
  desc = "Whether to print plots during analysis",
  option_name = "moo_print_plots",
  envvar_name = "MOO_PRINT_PLOTS"
)
options::define_option(
  option = "save_plots",
  default = TRUE,
  desc = "Whether to save plots to files during analysis",
  option_name = "moo_save_plots",
  envvar_name = "MOO_SAVE_PLOTS"
)
options::define_option(
  option = "plots_dir",
  default = "figures/",
  desc = "Path where plots are saved when `moo_save_plots` is `TRUE`",
  option_name = "moo_plots_dir",
  envvar_name = "MOO_PLOTS_DIR"
)


#' @eval options::as_roxygen_docs()
NULL

#' @title Option parameters
#' @eval options::as_params()
#' @name option_params
#'
NULL
