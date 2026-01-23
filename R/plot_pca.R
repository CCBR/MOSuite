#' Perform and plot a Principal Components Analysis
#'
#' @param moo_counts counts dataframe or `multiOmicDataSet` containing `count_type` & `sub_count_type` in the counts
#'   slot
#' @param principal_components vector with numbered principal components to plot. Use 2 for a 2D pca with ggplot, or 3
#'   for a 3D pca with plotly. (Default: `c(1,2)`)
#' @param ... additional arguments forwarded to method (see Details below)
#'
#' @examples
#' # multiOmicDataSet
#' moo <- multiOmicDataSet(
#'   sample_metadata = nidap_sample_metadata,
#'   anno_dat = data.frame(),
#'   counts_lst = list(
#'     "raw" = nidap_raw_counts,
#'     "clean" = nidap_clean_raw_counts
#'   )
#' )
#' plot_pca(moo, count_type = "clean", principal_components = c(1, 2))
#'
#' # 3D
#' plot_pca(moo, count_type = "clean", principal_components = c(1, 2, 3))
#'
#' # dataframe
#' plot_pca(nidap_clean_raw_counts,
#'   sample_metadata = nidap_sample_metadata,
#'   principal_components = c(1, 2)
#' )
#'
#' @details
#'
#'  See the low-level function docs for additional arguments
#'  depending on whether you're plotting 2 or 3 PCs:
#'
#'  - [plot_pca_2d] - used when there are **2** principal components
#'  - [plot_pca_3d] - used when there are **3** principal components
#'
#' # Methods
#'
#' | link to docs  | class  |
#' |---|---|
#' | [plot_pca_moo] | `multiOmicDataSet` |
#' | [plot_pca_dat] | `data.frame`       |
#'
#' @export
#' @return PCA plot (2D or 3D depending on the number of `principal_components`)
#'
#' @family plotters
#' @family PCA functions
#' @keywords plotters
#' @family moo methods
plot_pca <- S7::new_generic(
  "plot_pca",
  "moo_counts",
  function(moo_counts, principal_components = c(1, 2), ...) {
    return(S7::S7_dispatch())
  }
)

#' Plot 2D or 3D PCA for multiOmicDataset
#'
#' @param moo_counts `multiOmicDataSet` containing `count_type` & `sub_count_type` in the counts slot
#' @param count_type the type of counts to use. Must be a name in the counts slot (`names(moo@counts)`).
#' @param sub_count_type used if `count_type` is a list in the counts slot: specify the sub count type within the list.
#'   Must be a name in `names(moo@counts[[count_type]])`.
#' @param principal_components vector with numbered principal components to plot. Use 2 for a 2D pca with ggplot, or 3
#'   for a 3D pca with plotly. (Default: `c(1,2)`)
#' @param ... additional arguments forwarded to [plot_pca_2d()] (if 2 PCs) or [plot_pca_3d()] (if 3 PCs).
#'
#' @returns PCA plot
#'
#' @name plot_pca_moo
#' @seealso [plot_pca] generic
#' @family plotters for multiOmicDataSets
S7::method(plot_pca, multiOmicDataSet) <- function(
  moo_counts,
  count_type,
  sub_count_type = NULL,
  principal_components = c(1, 2),
  ...
) {
  counts_dat <- extract_counts(moo_counts, count_type, sub_count_type)
  return(
    plot_pca(
      counts_dat,
      sample_metadata = moo_counts@sample_meta,
      principal_components = principal_components,
      ...
    )
  )
}

#' Plot 2D or 3D PCA for counts dataframe
#'
#' @param moo_counts counts dataframe
#' @param sample_metadata **Required** if `moo_counts` is a `data.frame`: sample metadata as a data frame or tibble.
#' @param principal_components vector with numbered principal components to plot. Use 2 for a 2D pca with ggplot, or 3
#'   for a 3D pca with plotly. (Default: `c(1,2)`)
#' @param ... additional arguments forwarded to [plot_pca_2d()] (if 2 PCs) or [plot_pca_3d()] (if 3 PCs).
#'
#' @name plot_pca_dat
#' @seealso [plot_pca] generic
#' @family plotters for counts dataframes
S7::method(plot_pca, S7::class_data.frame) <- function(
  moo_counts,
  sample_metadata,
  principal_components = c(1, 2),
  ...
) {
  len_pcs <- length(principal_components)
  if (len_pcs == 2) {
    plot_fun <- plot_pca_2d
  } else if (len_pcs == 3) {
    plot_fun <- plot_pca_3d
  } else {
    stop(glue::glue(
      "Principal components must have exactly 2 or 3 items. Length: {len_pcs}"
    ))
  }
  return(
    plot_fun(
      moo_counts,
      sample_metadata = sample_metadata,
      principal_components = principal_components,
      ...
    )
  )
}

#' @rdname plot_pca_2d
#' @name plot_pca_2d
#' @export
plot_pca_2d <- S7::new_generic(
  "plot_pca_2d",
  "moo_counts",
  function(
    moo_counts,
    count_type = NULL,
    sub_count_type = NULL,
    sample_metadata = NULL,
    sample_id_colname = NULL,
    feature_id_colname = NULL,
    group_colname = "Group",
    label_colname = "Label",
    samples_to_rename = NULL,
    color_values = c(
      "#5954d6",
      "#e1562c",
      "#b80058",
      "#00c6f8",
      "#d163e6",
      "#00a76c",
      "#ff9287",
      "#008cf9",
      "#006e00",
      "#796880",
      "#FFA500",
      "#878500"
    ),
    principal_components = c(1, 2),
    legend_position = "top",
    point_size = 1,
    add_label = TRUE,
    label_font_size = 3,
    label_offset_x_ = 2,
    label_offset_y_ = 2,
    interactive_plots = FALSE,
    plots_subdir = "pca",
    plot_filename = "pca_2D.png",
    print_plots = options::opt("print_plots"),
    save_plots = options::opt("save_plots")
  ) {
    return(S7::S7_dispatch())
  }
)

#' @rdname plot_pca_2d
#' @name plot_pca_2d
#' @export
S7::method(plot_pca_2d, multiOmicDataSet) <- function(
  moo_counts,
  count_type = NULL,
  sub_count_type = NULL,
  sample_metadata = NULL,
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  group_colname = "Group",
  label_colname = "Label",
  samples_to_rename = NULL,
  color_values = c(
    "#5954d6",
    "#e1562c",
    "#b80058",
    "#00c6f8",
    "#d163e6",
    "#00a76c",
    "#ff9287",
    "#008cf9",
    "#006e00",
    "#796880",
    "#FFA500",
    "#878500"
  ),
  principal_components = c(1, 2),
  legend_position = "top",
  point_size = 1,
  add_label = TRUE,
  label_font_size = 3,
  label_offset_x_ = 2,
  label_offset_y_ = 2,
  interactive_plots = FALSE,
  plots_subdir = "pca",
  plot_filename = "pca_2D.png",
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots")
) {
  counts_dat <- extract_counts(moo_counts, count_type, sub_count_type)
  return(plot_pca_2d(
    counts_dat,
    sample_metadata = moo_counts@sample_meta,
    sample_id_colname = sample_id_colname,
    feature_id_colname = feature_id_colname,
    group_colname = group_colname,
    label_colname = label_colname,
    samples_to_rename = samples_to_rename,
    color_values = color_values,
    principal_components = principal_components,
    legend_position = legend_position,
    point_size = point_size,
    add_label = add_label,
    label_font_size = label_font_size,
    label_offset_x_ = label_offset_x_,
    label_offset_y_ = label_offset_y_,
    interactive_plots = interactive_plots,
    plots_subdir = plots_subdir,
    plot_filename = plot_filename,
    print_plots = print_plots,
    save_plots = save_plots
  ))
}

#' Perform and plot a 2D Principal Components Analysis
#'
#' @inheritParams create_multiOmicDataSet_from_dataframes
#' @inheritParams plot_histogram
#' @inheritParams plot_expr_heatmap
#' @inheritParams filter_counts
#'
#' @param sample_metadata sample metadata as a data frame or tibble.
#' @param sample_id_colname The column from the sample metadata containing the sample names. The names in this column
#'   must exactly match the names used as the sample column names of your input Counts Matrix. (Default: `NULL` - first
#'   column in the sample metadata will be used.)
#' @param feature_id_colname The column from the counts dataa containing the Feature IDs (Usually Gene or Protein ID).
#'   This is usually the first column of your input Counts Matrix. Only columns of Text type from your input Counts
#'   Matrix will be available to select for this parameter. (Default: `NULL` - first column in the counts matrix will be
#'   used.)
#' @param group_colname The column from the sample metadata containing the sample group information. This is usually a
#'   column showing to which experimental treatments each sample belongs (e.g. WildType, Knockout, Tumor, Normal,
#'   Before, After, etc.).
#' @param label_colname The column from the sample metadata containing the sample labels as you wish them to appear in
#'   the plots produced by this template. This can be the same Sample Names Column. However, you may desire different
#'   labels to display on your figure (e.g. shorter labels are sometimes preferred on plots). In that case, select the
#'   column with your preferred Labels here. The selected column should contain unique names for each sample. (Default:
#'   `NULL` -- `sample_id_colname` will be used.)
#' @param samples_to_rename If you do not have a Plot Labels Column in your sample metadata table, you can use this
#'   parameter to rename samples manually for display on the PCA plot. Use "Add item" to add each additional sample for
#'   renaming. Use the following format to describe which old name (in your sample metadata table) you want to rename to
#'   which new name: old_name: new_name
#' @param color_values vector of colors as hex values or names recognized by R
#' @param principal_components vector with numbered principal components to plot
#' @param legend_position passed to in `legend.position` `ggplot2::theme()`
#' @param point_size size for `ggplot2::geom_point()`
#' @param add_label whether to add text labels for the points
#'
#' @export
#' @return ggplot object
#'
#' @seealso [plot_pca] generic
#' @family PCA functions
#'
#' @rdname plot_pca_2d
#' @name plot_pca_2d
S7::method(plot_pca_2d, S7::class_data.frame) <- function(
  moo_counts,
  count_type = NULL,
  sub_count_type = NULL,
  sample_metadata = NULL,
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  group_colname = "Group",
  label_colname = "Label",
  samples_to_rename = NULL,
  color_values = c(
    "#5954d6",
    "#e1562c",
    "#b80058",
    "#00c6f8",
    "#d163e6",
    "#00a76c",
    "#ff9287",
    "#008cf9",
    "#006e00",
    "#796880",
    "#FFA500",
    "#878500"
  ),
  principal_components = c(1, 2),
  legend_position = "top",
  point_size = 1,
  add_label = TRUE,
  label_font_size = 3,
  label_offset_x_ = 2,
  label_offset_y_ = 2,
  interactive_plots = FALSE,
  plots_subdir = "pca",
  plot_filename = "pca_2D.png",
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots")
) {
  PC <- std.dev <- percent <- cumulative <- NULL
  if (length(principal_components) != 2) {
    stop(
      glue::glue(
        "principal_components must contain 2 values: {principal_components}"
      )
    )
  }

  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(moo_counts)[1]
  }

  # calculate PCA
  pca_df <- calc_pca(
    counts_dat = moo_counts,
    sample_metadata = sample_metadata,
    sample_id_colname = sample_id_colname,
    feature_id_colname = feature_id_colname
  ) %>%
    dplyr::filter(PC %in% principal_components) %>%
    # TODO consider redesigning to make rename_samples() unnecessary. Use Label column instead?
    rename_samples(samples_to_rename_manually = samples_to_rename)

  pca_wide <- pca_df %>%
    dplyr::select(-c(std.dev, percent, cumulative)) %>%
    tidyr::pivot_wider(
      names_from = "PC",
      names_prefix = "PC",
      values_from = "value"
    )
  prin_comp_x <- principal_components[1]
  prin_comp_y <- principal_components[2]
  # plot PCA
  pca_plot <- pca_wide %>%
    dplyr::mutate(
      !!rlang::sym(group_colname) := as.character(!!rlang::sym(group_colname))
    ) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = !!rlang::sym(glue::glue("PC{prin_comp_x}")),
      y = !!rlang::sym(glue::glue("PC{prin_comp_y}")),
      text = !!rlang::sym(sample_id_colname)
    )) +
    ggplot2::geom_point(
      ggplot2::aes(color = !!rlang::sym(group_colname)),
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
        show.legend = FALSE,
        direction = c("both"),
        box.padding = 1.25
      )
  }
  if (isTRUE(interactive_plots)) {
    pca_plot <- (pca_plot) %>%
      plotly::ggplotly(tooltip = c(sample_id_colname, group_colname))
  }

  print_or_save_plot(
    pca_plot,
    filename = file.path(plots_subdir, plot_filename),
    print_plots = print_plots,
    save_plots = save_plots
  )

  return(pca_plot)
}

#' @rdname plot_pca_3d
#' @name plot_pca_3d
#' @export
plot_pca_3d <- S7::new_generic(
  "plot_pca_3d",
  "moo_counts",
  function(
    moo_counts,
    count_type = NULL,
    sub_count_type = NULL,
    sample_metadata = NULL,
    feature_id_colname = NULL,
    sample_id_colname = NULL,
    samples_to_rename = NULL,
    group_colname = "Group",
    label_colname = "Label",
    principal_components = c(1, 2, 3),
    point_size = 8,
    label_font_size = 24,
    color_values = c(
      "#5954d6",
      "#e1562c",
      "#b80058",
      "#00c6f8",
      "#d163e6",
      "#00a76c",
      "#ff9287",
      "#008cf9",
      "#006e00",
      "#796880",
      "#FFA500",
      "#878500"
    ),
    plot_title = "PCA 3D",
    plot_filename = "pca_3D.html",
    print_plots = options::opt("print_plots"),
    save_plots = options::opt("save_plots"),
    plots_subdir = "pca"
  ) {
    return(S7::S7_dispatch())
  }
)

#' @rdname plot_pca_3d
#' @name plot_pca_3d
#' @export
S7::method(plot_pca_3d, multiOmicDataSet) <- function(
  moo_counts,
  count_type = NULL,
  sub_count_type = NULL,
  sample_metadata = NULL,
  feature_id_colname = NULL,
  sample_id_colname = NULL,
  samples_to_rename = NULL,
  group_colname = "Group",
  label_colname = "Label",
  principal_components = c(1, 2, 3),
  point_size = 8,
  label_font_size = 24,
  color_values = c(
    "#5954d6",
    "#e1562c",
    "#b80058",
    "#00c6f8",
    "#d163e6",
    "#00a76c",
    "#ff9287",
    "#008cf9",
    "#006e00",
    "#796880",
    "#FFA500",
    "#878500"
  ),
  plot_title = "PCA 3D",
  plot_filename = "pca_3D.html",
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = "pca"
) {
  counts_dat <- extract_counts(moo_counts, count_type, sub_count_type)
  return(
    plot_pca_3d(
      counts_dat,
      sample_metadata = moo_counts@sample_meta,
      count_type = count_type,
      sub_count_type = sub_count_type,
      feature_id_colname = feature_id_colname,
      sample_id_colname = sample_id_colname,
      samples_to_rename = samples_to_rename,
      group_colname = group_colname,
      label_colname = label_colname,
      principal_components = principal_components,
      point_size = point_size,
      label_font_size = label_font_size,
      color_values = color_values,
      plot_title = plot_title,
      plot_filename = plot_filename,
      print_plots = print_plots,
      save_plots = save_plots,
      plots_subdir = plots_subdir
    )
  )
}

#' 3D PCA for counts dataframe
#'
#' @inheritParams plot_pca_2d
#' @inheritParams filter_counts
#' @inheritParams plot_histogram
#' @inheritParams plot_expr_heatmap
#'
#' @param principal_components vector with numbered principal components to plot
#' @param point_size size for `ggplot2::geom_point()`
#' @param plot_title title for the plot
#'
#' @export
#' @returns `plotly::plot_ly` figure
#'
#' @family PCA functions
#'
#' @rdname plot_pca_3d
#' @name plot_pca_3d
S7::method(plot_pca_3d, S7::class_data.frame) <- function(
  moo_counts,
  count_type = NULL,
  sub_count_type = NULL,
  sample_metadata = NULL,
  feature_id_colname = NULL,
  sample_id_colname = NULL,
  samples_to_rename = NULL,
  group_colname = "Group",
  label_colname = "Label",
  principal_components = c(1, 2, 3),
  point_size = 8,
  label_font_size = 24,
  color_values = c(
    "#5954d6",
    "#e1562c",
    "#b80058",
    "#00c6f8",
    "#d163e6",
    "#00a76c",
    "#ff9287",
    "#008cf9",
    "#006e00",
    "#796880",
    "#FFA500",
    "#878500"
  ),
  plot_title = "PCA 3D",
  plot_filename = "pca_3D.html",
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = "pca"
) {
  PC <- std.dev <- percent <- cumulative <- NULL
  if (length(principal_components) != 3) {
    stop(
      glue::glue(
        "principal_components must contain 3 values: {principal_components}"
      )
    )
  }

  if (is.null(sample_metadata)) {
    stop("sample_metadata cannot be NULL")
  }
  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }

  # if (is.null(color_values)) {
  #   color_values <- moo_nidap@analyses[['colors']][['Group']]
  # }

  # calculate PCA
  pca_df <- calc_pca(
    counts_dat = moo_counts,
    sample_metadata = sample_metadata,
    sample_id_colname = sample_id_colname,
    feature_id_colname = feature_id_colname
  ) %>%
    dplyr::filter(PC %in% principal_components)
  pca_wide <- pca_df %>%
    dplyr::select(-c(std.dev, percent, cumulative)) %>%
    tidyr::pivot_wider(
      names_from = "PC",
      names_prefix = "PC",
      values_from = "value"
    )
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

  print_or_save_plot(
    fig,
    filename = file.path(plots_subdir, plot_filename),
    print_plots = print_plots,
    save_plots = save_plots
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

#' Perform principal components analysis
#'
#' @param counts_dat data frame of feature counts (e.g. from the counts slot of a `multiOmicDataSet`).
#' @param sample_metadata sample metadata as a data frame or tibble.
#' @param sample_id_colname The column from the sample metadata containing the sample names. The names in this column
#'   must exactly match the names used as the sample column names of your input Counts Matrix. (Default: `NULL` - first
#'   column in the sample metadata will be used.)
#' @param feature_id_colname The column from the counts dataa containing the Feature IDs (Usually Gene or Protein ID).
#'   This is usually the first column of your input Counts Matrix. Only columns of Text type from your input Counts
#'   Matrix will be available to select for this parameter. (Default: `NULL` - first column in the counts matrix will be
#'   used.)
#'
#' @returns data frame with statistics for each principal component
#' @export
#'
#' @examples
#' calc_pca(nidap_raw_counts, nidap_sample_metadata) %>% head()
#' @family PCA functions
calc_pca <- function(
  counts_dat,
  sample_metadata,
  sample_id_colname = NULL,
  feature_id_colname = NULL
) {
  var <- row <- percent <- NULL
  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(counts_dat)[1]
  }
  counts_dat <- counts_dat %>%
    as.data.frame() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(feature_id_colname)
  # sample-wise PCA
  tedf <- t(counts_dat)
  # remove samples with all NAs
  tedf_filt <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
  # remove samples with zero variance
  tedf_var <- tedf_filt[, apply(tedf_filt, 2, var) != 0]
  # calculate PCA
  pca_fit <- stats::prcomp(tedf_var, scale = TRUE)
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
