calc_pca <- function(counts_dat,
                     sample_metadata,
                     sample_id_colname = NULL) {
  var <- xdata <- ydata <- group <- row <- NULL
  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (sample_id_colname %in% counts_dat) {
    counts_dat %<>% tibble::column_to_rownames(sample_id_colname)
  }
  # sample-wise PCA
  tedf <- t(counts_dat)
  # remove samples with all NAs
  tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
  # remove samples with zero variance
  tedf_var <- tedf[, apply(tedf, 2, var) != 0]
  # calculate PCA
  pca_fit <- stats::prcomp(tedf_var, scale = T)
  pca_df <- pca_fit %>%
    broom::tidy() %>%
    dplyr::rename(!!rlang::sym(sample_id_colname) := row) %>%
    dplyr::left_join(pca_fit %>%
      broom::tidy(matrix = "eigenvalues") %>%
      dplyr::mutate(percent = percent * 100)) %>%
    dplyr::left_join(sample_metadata)
  return(pca_df)
}

#' Perform and plot a Principal Components Analysis
#'
#' @inheritParams filter_counts
#' @inheritParams plot_histogram
#'
#' @return ggplot object
#'
plot_pca <- function(counts_dat,
                     sample_metadata,
                     sample_id_colname = NULL,
                     samples_to_rename = NULL,
                     group_colname = "Group",
                     label_colname = "Label",
                     color_values = c(),
                     principal_components = c(1, 2),
                     legend_position_for_pca = "top",
                     point_size_for_pca = 1,
                     add_label_to_pca = TRUE,
                     label_font_size = 3,
                     label_offset_x_ = 2,
                     label_offset_y_ = 2,
                     make_plots_interactive = FALSE) {
  if (length(principal_components) < 2 || length(principal_components) > 3) {
    stop(glue::glue("principal_components must contain 2 or 3 values: {principal_components}"))
  }

  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }

  # calculate PCA
  pca_df <- calc_pca(
    counts_dat = counts_dat,
    sample_metadata = sample_metadata,
    sample_id_colname = sample_id_colname
  ) %>%
    dplyr::filter(PC %in% principal_components) %>%
    # TODO consider redesigning to make this unnecessary. Use Label column instead?
    rename_samples(samples_to_rename_manually = samples_to_rename)
  pca_wide <- pca_df %>%
    select(-c(std.dev, percent, cumulative)) %>%
    tidyr::pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value")
  prin_comp_x <- principal_components[1]
  prin_comp_y <- principal_components[2]
  # plot PCA
  pca_plot <- pca_wide %>%
    ggplot2::ggplot(ggplot2::aes(
      x = !!rlang::sym(glue::glue("PC{prin_comp_x}")),
      y = !!rlang::sym(glue::glue("PC{prin_comp_y}")),
      text = !!rlang::sym(sample_id_colname)
    )) +
    ggplot2::geom_point(ggplot2::aes(color = !!rlang::sym(group_colname)),
      size = point_size_for_pca
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = legend_position_for_pca,
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

  if (add_label_to_pca == TRUE) {
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
  if (isTRUE(make_plots_interactive)) {
    pca_plot <- (pca_plot) %>%
      plotly::ggplotly(tooltip = c(sample_id_colname, group_colname))
  }
  return(pca_plot)
}

get_pc_percent_lab <- function(pca_df, pc) {
  PC <- NULL
  perc <- pca_df %>%
    dplyr::filter(PC == pc) %>%
    dplyr::pull(percent) %>%
    unique() %>%
    round(digits = 1)
  return(glue::glue("PC{pc} {perc}%"))
}
