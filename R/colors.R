utils::globalVariables("mosuite_palette")

#' Get random colors.
#'
#' Note: this function is not guaranteed to create a color blind friendly
#' palette. Consider using other palettes such as
#' `RColorBrewer::display.brewer.all(colorblindFriendly = TRUE)`.
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
#' Returns non-missing values from `dat[[colname]]`. For factor columns, values
#' are returned in factor-level order; otherwise, values keep first-observed
#' order.
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
#' @param palette Character vector of colors to assign. Defaults to
#'   `mosuite_palette`.
#'
#' @returns named list, with each column in `sample_metadata` containing a corresponding entry with a named vector of
#'   colors
#' @export
#'
#' @examples
#' get_colors_lst(nidap_sample_metadata)
#' get_colors_lst(nidap_sample_metadata, palette = RColorBrewer::brewer.pal(12, "Set3"))
get_colors_lst <- function(
  sample_metadata,
  palette = mosuite_palette
) {
  dat_colnames <- colnames(sample_metadata)
  n_palette <- length(palette)

  color_offset <- 0L
  color_lists <- vector("list", length(dat_colnames))
  for (i in seq_along(dat_colnames)) {
    colname <- dat_colnames[[i]]
    n_obs <- length(get_observed_values(sample_metadata, colname))
    # Only offset when the column is small enough that unique colors are available
    use_offset <- n_obs <= n_palette / 2 && color_offset + n_obs <= n_palette
    vctr <- get_colors_vctr(
      dat = sample_metadata,
      colname = colname,
      palette = palette,
      color_offset = if (use_offset) color_offset else 0L
    )
    if (use_offset) {
      color_offset <- color_offset + n_obs
    }
    color_lists[[i]] <- vctr
  }
  names(color_lists) <- dat_colnames
  return(color_lists)
}

#' Get vector of colors for observations in one column of a data frame
#'
#' Assigns one color per unique observed value in `dat[[colname]]`, drawn from
#' `palette` starting at `color_offset`. If the palette is too short,
#' falls back to `get_random_colors()`. Factor columns use factor-level order;
#' other columns use first-observed order.
#'
#' @inheritParams get_colors_lst
#' @param dat data frame
#' @param colname column name in `dat`
#' @param color_offset integer; number of palette colors to skip before
#'   assigning colors to this column's values. Used by [get_colors_lst()] to
#'   avoid repeating colors across columns with few unique values.
#' @returns Named character vector of hex colors, one per unique observed value
#'   in `dat[[colname]]`.
#' @export
#'
#' @examples
#' get_colors_vctr(nidap_sample_metadata, "Group")
#' get_colors_vctr(nidap_sample_metadata, "Group", color_offset = 3L)
get_colors_vctr <- function(
  dat,
  colname,
  palette = mosuite_palette,
  color_offset = 0L
) {
  obs <- get_observed_values(dat, colname)
  n_obs <- length(obs)

  if (n_obs == 0) {
    colors_vctr <- character(0)
  } else {
    # if fewer colors are available than needed, handle gracefully
    if (length(palette) < n_obs + color_offset) {
      # If an offset pushed us past the palette end but n_obs alone would fit,
      # retry from the start of the palette before falling back to random colors.
      if (color_offset > 0L && length(palette) >= n_obs) {
        color_offset <- 0L
      }
    }
    # If still not enough colors, fall back to random
    if (length(palette) < n_obs) {
      message(glue::glue(
        "Number of unique values ({n_obs}) in column \"{colname}\" ",
        "exceeds the palette maximum. Falling back to random colors."
      ))
      colors_vctr <- get_random_colors(n_obs)
    } else {
      colors_vctr <- palette[seq.int(color_offset + 1L, color_offset + n_obs)]
    }
    names(colors_vctr) <- obs
  }

  return(colors_vctr)
}


#' Resolve plotting colors for one column
#'
#' Uses `color_values` when supplied; otherwise generates colors with
#' [get_colors_vctr()]. If `color_values` is named and covers all observed
#' values, it is returned as-is. If too few colors are provided, missing colors
#' are generated and appended.
#'
#' @param dat data frame
#' @param colname column name in `dat`
#' @param color_values optional named or unnamed character vector of colors
#' @param palette character vector of colors used to generate defaults
#' @returns named character vector of colors matching observed values in
#'   `dat[[colname]]`
#' @keywords internal
resolve_plot_colors <- function(
  dat,
  colname,
  color_values = NULL,
  palette = mosuite_palette
) {
  obs <- get_observed_values(dat, colname)

  if (length(obs) == 0) {
    return(color_values)
  }

  if (is.null(color_values)) {
    return(get_colors_vctr(dat, colname, palette = palette))
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
      palette = palette
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

  p <- plot_palette(df) +
    ggplot2::labs(title = "mosuite_palette")

  return(p)
}

#' Display colors for a multiOmicDataSet object
#'
#' Plots a palette strip for each group column stored in `moo@analyses$colors`,
#' stacked vertically. Each strip shows the assigned hex colors and their codes.
#'
#' @param moo A `multiOmicDataSet` object (see
#'   [create_multiOmicDataSet_from_dataframes()]).
#' @returns A [patchwork][patchwork::wrap_plots] of [ggplot2::ggplot] objects,
#'   one per group column in `moo@analyses$colors`.
#' @export
#' @examples
#' moo <- create_multiOmicDataSet_from_dataframes(nidap_sample_metadata, nidap_raw_counts)
#' display_colors(moo)
display_colors <- function(moo) {
  abort_packages_not_installed("patchwork")
  colors_lst <- moo@analyses$colors
  plots_lst <- lapply(names(colors_lst), function(colname) {
    palette <- colors_lst[[colname]]
    n <- length(palette)
    df <- data.frame(
      hex = palette,
      idx = factor(seq_len(n))
    )
    return(
      p <- plot_palette(df) +
        ggplot2::labs(title = colname)
    )
  })
  return(patchwork::wrap_plots(plots_lst, ncol = 1))
}

#' Plot a palette tile strip
#'
#' Renders a data frame with columns `hex` and `idx` as a row of colored tiles,
#' each labeled with its hex code. Used internally by [display_palette()] and
#' [display_colors()].
#'
#' @param dat data frame with columns `hex` (hex color codes) and `idx`
#'   (factor, used for faceting)
#' @returns a [ggplot2::ggplot] object
#' @keywords internal
plot_palette <- function(dat) {
  p <- ggplot2::ggplot(dat) +
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
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
  return(p)
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
#' @returns `moo` with colors updated at `moo@analyses$colors[[colname]]`
#' @export
#'
#' @examples
#' moo <- create_multiOmicDataSet_from_dataframes(
#'   sample_metadata = as.data.frame(nidap_sample_metadata),
#'   counts_dat = as.data.frame(nidap_raw_counts)
#' )
#' moo@analyses$colors$Group
#' moo <- moo |> set_color_pal("Group", palette = RColorBrewer::brewer.pal(3, "Set2"))
#' moo@analyses$colors$Group
#'
#' @family moo methods
set_color_pal <- S7::new_generic(
  "set_color_pal",
  "moo",
  function(moo, colname, palette = mosuite_palette) {
    return(S7::S7_dispatch())
  }
)

S7::method(set_color_pal, multiOmicDataSet) <- function(
  moo,
  colname,
  palette = mosuite_palette
) {
  moo@analyses$colors[[colname]] <- get_colors_vctr(
    dat = moo@sample_meta,
    colname = colname,
    palette = palette
  )
  return(moo)
}
