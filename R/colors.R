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
  ourColorSpace <- colorspace::RGB(
    stats::runif(n),
    stats::runif(n),
    stats::runif(n)
  )
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
#' @param palette_fun Function for selecting colors. Assumed to contain `n` for the number of colors. Default:
#'   `grDevices::palette.colors()`
#' @param ... additional arguments forwarded to `palette_fun`
#'
#' @returns named list, with each column in `sample_metadata` containing entry with a named vector of colors
#' @export
#'
#' @examples
#' get_colors_lst(nidap_sample_metadata)
#' \dontrun{
#' get_colors_lst(nidap_sample_metadata, palette_fun = RColorBrewer::brewer.pal, name = "Set3")
#' }
get_colors_lst <- function(
  sample_metadata,
  palette_fun = grDevices::palette.colors,
  ...
) {
  dat_colnames <- colnames(sample_metadata)
  color_lists <- dat_colnames |>
    purrr::map(
      .f = get_colors_vctr,
      dat = sample_metadata,
      palette_fun = palette_fun,
      ...
    )
  names(color_lists) <- dat_colnames
  return(color_lists)
}

#' Get vector of colors for observations in one column of a data frame
#'
#' @inheritParams get_colors_lst
#' @param dat data frame
#' @param colname column name in `dat`
#' @returns named vector of colors for each unique observation in `dat$colname`
#' @export
#'
get_colors_vctr <- function(
  dat,
  colname,
  palette_fun = grDevices::palette.colors,
  ...
) {
  obs <- dat |>
    dplyr::pull(colname) |>
    unique()
  n_obs <- length(obs)

  warned_cnd <- NULL
  colors_vctr <- withCallingHandlers(
    warning = function(cnd) {
      warned_cnd <<- cnd
      invokeRestart("muffleWarning")
    },
    palette_fun(n = n_obs, ...)
  )

  # if fewer colors were returned than needed (e.g. when n exceeds the palette maximum,
  # such as Okabe-Ito's maximum of 9), silently fall back to random colors
  if (length(colors_vctr) < n_obs) {
    colors_vctr <- get_random_colors(n_obs)
  } else if (!is.null(warned_cnd)) {
    # warning was raised but we still have enough colors (e.g. brewer.pal warns when n < 3
    # but returns 3 colors); convert to a message and re-raise the original warning
    message(glue::glue(
      'Warning raised in get_color_vctr() for column "{colname}"'
    ))
    warning(conditionMessage(warned_cnd))
  }

  # if more colors are returned than are in the observations, truncate the vector.
  # this occurs when using RColorBrewer::brewer.pal with n < 3
  colors_vctr <- colors_vctr[seq_len(n_obs)]

  names(colors_vctr) <- obs
  return(colors_vctr)
}

#' Set color palette for a single group/column
#'
#' This allows you to set custom palettes individually for groups in the dataset
#'
#' @inheritParams get_colors_lst
#'
#' @param moo `multiOmicDataSet` object (see `create_multiOmicDataSet_from_dataframes()`)
#' @param colname group column name to set the palette for
#'
#' @returns `moo` with colors updated at `moo@analyses$colors$colname`
#' @export
#'
#' @examples
#' moo <- create_multiOmicDataSet_from_dataframes(
#'   sample_metadata = as.data.frame(nidap_sample_metadata),
#'   counts_dat = as.data.frame(nidap_raw_counts)
#' )
#' moo@analyses$colors$Group
#' moo <- moo |> set_color_pal("Group", palette_fun = RColorBrewer::brewer.pal, name = "Set2")
#' moo@analyses$colors$Group
#'
#' @family moo methods
set_color_pal <- S7::new_generic(
  "set_color_pal",
  "moo",
  function(moo, colname, palette_fun = grDevices::palette.colors, ...) {
    return(S7::S7_dispatch())
  }
)

S7::method(set_color_pal, multiOmicDataSet) <- function(
  moo,
  colname,
  palette_fun = grDevices::palette.colors,
  ...
) {
  moo@analyses[["colors"]][[colname]] <- get_colors_vctr(
    dat = moo@sample_meta,
    colname = colname,
    palette_fun = palette_fun,
    ...
  )
  return(moo)
}
