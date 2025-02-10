#' Get random colors.
#'
#' Note: this function is not guaranteed to create a color blind friendly palette.
#' Consider using other palettes such as `RColorBrewer::display.brewer.all(colorblindFriendly = TRUE)`.
#'
#' @param num_colors number of colors to select.
#' @param n number of random RGB values to generate in the color space.
#'
#' @return vector of random colors in hex format.
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' set.seed(10)
#' get_random_colors(5)
#' }
get_random_colors <- function(num_colors, n = 2e3) {
  abort_packages_not_installed("colorspace")
  if (num_colors < 1) {
    stop("num_colors must be at least 1")
  }
  n <- 2e3
  ourColorSpace <- colorspace::RGB(stats::runif(n), stats::runif(n), stats::runif(n))
  ourColorSpace <- methods::as(ourColorSpace, "LAB")
  currentColorSpace <- ourColorSpace@coords
  # Set iter.max to 20 to avoid convergence warnings.
  km <- stats::kmeans(currentColorSpace, num_colors, iter.max = 20)
  return(unname(colorspace::hex(colorspace::LAB(km$centers))))
}


#' Create named list of default colors for plotting
#'
#' @inheritParams create_multiOmicDataSet_from_dataframes
#'
#' @param palette_fun Function for selecting colors. Assumed to contain `n` for the number of colors. Default: `grDevices::palette.colors()`
#' @param ... additional arguments forwarded to `palette_fun`
#'
#' @returns named list, with each column in `sample_metadata` containing entry with a named vector of colors
#' @export
#'
#' @examples
#' set_colors(nidap_sample_metadata)
#' \dontrun{
#' set_colors(nidap_sample_metadata, palette_fun = RColorBrewer::brewer.pal, name = "Set3")
#' }
set_colors <- function(sample_metadata, palette_fun = grDevices::palette.colors, ...) {
  dat_colnames <- colnames(sample_metadata)
  color_lists <- dat_colnames %>% lapply(function(colname, dat = sample_metadata) {
    obs <- dat %>%
      dplyr::pull(colname) %>%
      unique()
    withCallingHandlers(
      warning = function(cnd) {
        message(glue::glue('Warning raised in set_colors() for column "{colname}"'))
      },
      colors_vctr <- palette_fun(n = length(obs), ...)
    )

    # if more colors are returned than are in the observations, truncate the vector.
    # this occurs when using RColorBrewer::brewer.pal with n < 3
    colors_vctr <- colors_vctr[1:length(obs)]

    names(colors_vctr) <- obs
    return(colors_vctr)
  })
  names(color_lists) <- dat_colnames
  return(color_lists)
}
