utils::globalVariables("mosuite_palette")

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

#' Select colors from MOSuite's default palette
#'
#' @param n number of colors to select.
#' @param ... additional arguments (ignored).
#'
#' @returns vector of colors in hex format.
#' @export
#'
#' @examples
#' select_mosuite_colors(5)
select_mosuite_colors <- function(n, ...) {
  return(mosuite_palette[seq_len(min(n, length(mosuite_palette)))])
}

#' Get observed values from a column
#'
#' Returns non-missing values from `dat[[colname]]`. For factor columns, values are
#' returned in factor-level order; otherwise, values keep first-observed order.
#'
#' @param dat data frame
#' @param colname column name in `dat`
#' @returns character vector of observed values
#' @keywords internal
get_observed_values <- function(dat, colname) {
  values <- dplyr::pull(dat, colname)
  observed_values <- stats::na.omit(as.character(values))

  if (is.factor(values)) {
    return(levels(values)[levels(values) %in% observed_values])
  }

  return(unique(observed_values))
}


#' Create named list of default colors for plotting
#'
#' @inheritParams create_multiOmicDataSet_from_dataframes
#'
#' @param palette_fun Function for selecting colors. Assumed to contain `n` for the number of colors. Defaults to
#'   MOSuite's default plot palette. To use the previous R default palette behavior, pass
#'   `grDevices::palette.colors`.
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
  palette_fun = select_mosuite_colors,
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
#' @returns named vector of colors for each unique observation in `dat$colname`. Factor columns use factor level order;
#'   other columns use first-observed order.
#' @export
#'
get_colors_vctr <- function(
  dat,
  colname,
  palette_fun = select_mosuite_colors,
  ...
) {
  obs <- get_observed_values(dat, colname)
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
  # such as Okabe-Ito's maximum of 9), fall back to random colors
  if (length(colors_vctr) < n_obs) {
    message(glue::glue(
      "Number of unique values ({n_obs}) in column \"{colname}\" ",
      "exceeds the palette maximum. Falling back to random colors."
    ))
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


#' Resolve plotting colors for one column
#'
#' Uses `color_values` when supplied; otherwise generates colors with
#' `get_colors_vctr()`. If too few colors are provided, missing colors are generated
#' and appended.
#'
#' @param dat data frame
#' @param colname column name in `dat`
#' @param color_values optional vector of colors
#' @param palette_fun function used to generate colors
#' @param ... additional arguments forwarded to `palette_fun`
#' @returns named vector of colors matching observed values in `colname`
#' @keywords internal
resolve_plot_colors <- function(
  dat,
  colname,
  color_values = NULL,
  palette_fun = select_mosuite_colors,
  ...
) {
  obs <- get_observed_values(dat, colname)

  if (length(obs) == 0) {
    return(color_values)
  }

  if (is.null(color_values)) {
    return(get_colors_vctr(dat, colname, palette_fun = palette_fun, ...))
  }

  if (!is.null(names(color_values))) {
    if (all(obs %in% names(color_values))) {
      return(color_values)
    }
  }

  if (length(color_values) < length(obs)) {
    n_missing <- length(obs) - length(color_values)
    message(glue::glue(
      "color_values contains {length(color_values)} colors for ",
      "{length(obs)} values in column {colname}. Generating ",
      "{n_missing} additional colors."
    ))
    generated_colors <- get_colors_vctr(
      dat,
      colname,
      palette_fun = palette_fun,
      ...
    )
    color_values <- c(
      unname(color_values),
      unname(generated_colors)[seq.int(length(color_values) + 1, length(obs))]
    )
  }

  return(stats::setNames(unname(color_values)[seq_along(obs)], obs))
}

#' Display the mosuite color palette
#'
#' Plots each color in `mosuite_palette` as a labeled tile with its hex code
#' displayed below. The plot is rendered at a width proportional to the number
#' of colors so labels remain horizontal and legible.
#'
#' @param palette Character vector of hex color codes. Defaults to
#'   `mosuite_palette`.
#'
#' @returns Invisibly returns the underlying [ggplot2::ggplot] object.
#' @export
#'
#' @examples
#' display_palette()
#' display_palette(c("#FF0000", "#00FF00", "#0000FF"))
display_palette <- function(palette = mosuite_palette) {
  n <- length(palette)

  df <- data.frame(
    hex = palette,
    idx = factor(seq_len(n))
  )

  p <- ggplot2::ggplot(df) +
    ggplot2::geom_rect(
      ggplot2::aes(fill = .data$hex),
      xmin = 0,
      xmax = 1,
      ymin = 0.25,
      ymax = 1
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$hex),
      x = 0.5,
      y = 0.12,
      size = 2.8,
      family = "mono",
      vjust = 1
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::facet_wrap(~idx, nrow = 1) +
    ggplot2::theme_void() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(
        size = 9,
        margin = ggplot2::margin(b = 3)
      ),
      panel.spacing = ggplot2::unit(3, "pt"),
      plot.title = ggplot2::element_text(
        size = 12,
        margin = ggplot2::margin(b = 8)
      ),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    ) +
    ggplot2::labs(title = "mosuite_palette")

  # Render at a width that gives each tile enough room for its horizontal label
  panel_w_in <- 0.85
  total_w <- n * panel_w_in + 0.4
  total_h <- 2.5

  tmp <- tempfile(fileext = ".png")
  ggplot2::ggsave(tmp, plot = p, width = total_w, height = total_h, dpi = 150)

  img <- png::readPNG(tmp)
  grid::grid.newpage()
  grid::grid.draw(grid::rasterGrob(img, width = grid::unit(1, "npc")))

  return(invisible(p))
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
  function(moo, colname, palette_fun = select_mosuite_colors, ...) {
    return(S7::S7_dispatch())
  }
)

S7::method(set_color_pal, multiOmicDataSet) <- function(
  moo,
  colname,
  palette_fun = select_mosuite_colors,
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
