#' Perform principal components analysis
#'
#' @inheritParams plot_pca
#'
#' @returns data frame with statistics for each principal component
#' @export
#'
#' @examples
#' calc_pca(nidap_raw_counts, nidap_sample_metadata) %>% head()
calc_pca <- function(counts_dat,
                     sample_metadata,
                     sample_id_colname = NULL,
                     feature_id_colname = NULL) {
  var <- xdata <- ydata <- group <- row <- percent <- NULL
  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(counts_dat)[1]
  }
  counts_dat %<>% as.data.frame() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(feature_id_colname)
  # sample-wise PCA
  tedf <- t(counts_dat)
  # remove samples with all NAs
  tedf_filt <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
  # remove samples with zero variance
  tedf_var <- tedf_filt[, apply(tedf_filt, 2, var) != 0]
  # calculate PCA
  pca_fit <- stats::prcomp(tedf_var, scale = T)
  pca_df <- pca_fit %>%
    broom::tidy() %>%
    dplyr::rename(!!rlang::sym(sample_id_colname) := row) %>%
    dplyr::left_join(
      pca_fit %>%
        broom::tidy(matrix = "eigenvalues") %>%
        dplyr::mutate(percent = percent * 100),
      by = "PC"
    ) %>%
    dplyr::left_join(sample_metadata, by = sample_id_colname)
  return(pca_df)
}

#' Perform and plot a Principal Components Analysis
#'
#' @inheritParams create_multiOmicDataSet_from_dataframes
#' @inheritParams plot_histogram
#' @inheritParams filter_counts
#'
#' @param principal_components vector with numbered principal components to plot (Default: `c(1,2)`)
#' @param point_size size for `ggplot2::geom_point()`
#' @param add_label whether to add text labels for the points
#' @export
#'
#' @return ggplot object
#' @examples
#' plot_pca(nidap_raw_counts, nidap_sample_metadata)
plot_pca <- function(counts_dat,
                     sample_metadata,
                     sample_id_colname = NULL,
                     feature_id_colname = NULL,
                     samples_to_rename = NULL,
                     group_colname = "Group",
                     label_colname = "Label",
                     color_values = c(
                       "#5954d6", "#e1562c", "#b80058", "#00c6f8", "#d163e6", "#00a76c",
                       "#ff9287", "#008cf9", "#006e00", "#796880", "#FFA500", "#878500"
                     ),
                     principal_components = c(1, 2),
                     legend_position = "top",
                     point_size = 1,
                     add_label = TRUE,
                     label_font_size = 3,
                     label_offset_x_ = 2,
                     label_offset_y_ = 2,
                     interactive_plots = FALSE) {
  PC <- std.dev <- percent <- cumulative <- NULL
  if (length(principal_components) != 2) {
    stop(glue::glue("principal_components must contain 2 values: {principal_components}"))
  }

  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(counts_dat)[1]
  }

  # calculate PCA
  pca_df <- calc_pca(
    counts_dat = counts_dat,
    sample_metadata = sample_metadata,
    sample_id_colname = sample_id_colname,
    feature_id_colname = feature_id_colname
  ) %>% dplyr::filter(PC %in% principal_components) %>%
    # TODO consider redesigning to make rename_samples() unnecessary. Use Label column instead?
    rename_samples(samples_to_rename_manually = samples_to_rename)

  pca_wide <- pca_df %>%
    dplyr::select(-c(std.dev, percent, cumulative)) %>%
    tidyr::pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value")
  prin_comp_x <- principal_components[1]
  prin_comp_y <- principal_components[2]
  # plot PCA
  pca_plot <- pca_wide %>%
    dplyr::mutate(!!rlang::sym(group_colname) := as.character(!!rlang::sym(group_colname))) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = !!rlang::sym(glue::glue("PC{prin_comp_x}")),
      y = !!rlang::sym(glue::glue("PC{prin_comp_y}")),
      text = !!rlang::sym(sample_id_colname)
    )) +
    ggplot2::geom_point(ggplot2::aes(color = !!rlang::sym(group_colname)),
      size = point_size
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = legend_position,
      legend.title = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 18),
      axis.title = ggplot2::element_text(size = 20),
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1
      ),
      axis.ticks = ggplot2::element_line(linewidth = 1),
      legend.text = ggplot2::element_text(size = 18)
    ) +
    ggplot2::coord_fixed(ratio = 1.5) +
    ggplot2::scale_colour_manual(values = color_values) +
    ggplot2::xlab(get_pc_percent_lab(pca_df, prin_comp_x)) +
    ggplot2::ylab(get_pc_percent_lab(pca_df, prin_comp_y))

  if (add_label == TRUE) {
    pca_plot <- pca_plot +
      ggrepel::geom_text_repel(
        ggplot2::aes(
          label = !!rlang::sym(label_colname),
          color = !!rlang::sym(group_colname)
        ),
        size = 7,
        show.legend = F,
        direction = c("both"),
        box.padding = 1.25
      )
  }
  if (isTRUE(interactive_plots)) {
    pca_plot <- (pca_plot) %>%
      plotly::ggplotly(tooltip = c(sample_id_colname, group_colname))
  }
  return(pca_plot)
}

#' Plot 3-Dimensional PCA with plotly
#'
#' @inheritParams plot_pca
#' @inheritParams filter_counts
#'
#' @param principal_components vector with numbered principal components to plot (Default: `c(1,2,3)`)
#' @param plot_title title for the plot
#'
#' @returns `plotly::plot_ly` figure
#' @export
#'
#' @examples
#' plot_pca_3d(nidap_raw_counts, nidap_sample_metadata)
plot_pca_3d <- function(counts_dat,
                        sample_metadata,
                        sample_id_colname = NULL,
                        samples_to_rename = NULL,
                        group_colname = "Group",
                        label_colname = "Label",
                        principal_components = c(1, 2, 3),
                        point_size = 8,
                        label_font_size = 24,
                        color_values = c(
                          "#5954d6", "#e1562c", "#b80058", "#00c6f8", "#d163e6", "#00a76c",
                          "#ff9287", "#008cf9", "#006e00", "#796880", "#FFA500", "#878500"
                        ),
                        plot_title = "PCA 3D") {
  PC <- std.dev <- percent <- cumulative <- NULL
  if (length(principal_components) != 3) {
    stop(glue::glue("principal_components must contain 3 values: {principal_components}"))
  }

  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }

  # if (is.null(color_values)) {
  #   color_values <- moo_nidap@analyses[['colors']][['Group']]
  # }

  # calculate PCA
  pca_df <- calc_pca(
    counts_dat = counts_dat,
    sample_metadata = sample_metadata,
    sample_id_colname = sample_id_colname
  ) %>%
    dplyr::filter(PC %in% principal_components)
  pca_wide <- pca_df %>%
    dplyr::select(-c(std.dev, percent, cumulative)) %>%
    tidyr::pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value")
  prin_comp_x <- principal_components[1]
  prin_comp_y <- principal_components[2]
  prin_comp_z <- principal_components[3]

  fig <- plotly::plot_ly(
    pca_wide,
    x = stats::as.formula(paste0("~ PC", prin_comp_x)),
    y = stats::as.formula(paste0("~ PC", prin_comp_y)),
    z = stats::as.formula(paste0("~ PC", prin_comp_z)),
    color = stats::as.formula(paste("~", group_colname)),
    colors = color_values,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = point_size),
    hoverinfo = "text",
    text = stats::as.formula(paste("~", sample_id_colname)),
    size = label_font_size
  )
  return(fig)
}

#' Get label for Principal Component with percent of variation
#'
#' @param pca_df data frame from `calc_pca()`
#' @param pc which principal component to report (e.g. `1`)
#'
#' @returns glue string formatted with PC's percent of variation
#' @keywords internal
#' @examples
#' \dontrun{
#' data.frame(PC = c(1, 2, 3), percent = c(40, 10, 0.5)) %>%
#'   get_pc_percent_lab(2)
#' }
get_pc_percent_lab <- function(pca_df, pc) {
  PC <- percent <- NULL
  perc <- pca_df %>%
    dplyr::filter(PC == pc) %>%
    dplyr::pull(percent) %>%
    unique() %>%
    round(digits = 1)
  return(glue::glue("PC{pc} {perc}%"))
}
