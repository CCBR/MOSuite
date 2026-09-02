#' Filter low counts
#'
#' This is often the first step in the QC portion of an analysis to filter out
#' features that have very low raw counts across most or all of your samples.
#'
#' This function takes a multiOmicDataSet containing clean raw counts and a sample
#' metadata table, and returns the multiOmicDataSet object with filtered counts.
#' It also produces an image consisting of three QC plots.
#'
#' You can tune the threshold for tuning how low counts for a given gene are
#' before they are deemed "too low" and filtered out of downstream analysis. By
#' default, this parameter is set to 1, meaning any raw count value less than 1
#' will count as "too low".
#'
#' The QC plots are provided to help you assess: (1) PCA Plot: the within and
#' between group variance in expression after dimensionality reduction; (2)
#' Count Density Histogram: the dis/similarity of count distributions between
#' samples; and (3) Similarity Heatmap: the overall similarity of samples to one
#' another based on unsupervised clustering.
#'
#' @inheritParams option_params
#'
#' @param moo multiOmicDataSet object (see `create_multiOmicDataSet_from_dataframes()`)
#' @param count_type the type of counts to use -- must be a name in the counts slot (`moo@counts`)
#' @param feature_id_colname The column from the counts data containing the Feature IDs (Usually Gene or Protein ID).
#'   This is usually the first column of your input Counts Matrix. Only columns of Text type from your input Counts
#'   Matrix will be available to select for this parameter. (Default: `NULL` - first column in the counts matrix will be
#'   used.)
#' @param sample_id_colname The column from the sample metadata containing the sample names. The names in this column
#'   must exactly match the names used as the sample column names of your input Counts Matrix. (Default: `NULL` - first
#'   column in the sample metadata will be used.)
#' @param group_colname The column from the sample metadata containing the sample group information. This is usually a
#'   column showing to which experimental treatments each sample belongs (e.g. WildType, Knockout, Tumor, Normal,
#'   Before, After, etc.).
#' @param label_colname The column from the sample metadata containing the sample labels as you wish them to appear in
#'   heatmap and PCA figures. This can be the same Sample Names Column. However, you may desire different labels to
#'   display on your figures (e.g. shorter labels are sometimes preferred on plots). In that case, select the column with
#'   your preferred Labels here. The selected column should contain unique names for each sample. Use `add_label_to_pca`
#'   to control whether these labels are displayed on the PCA plot.
#' @param samples_to_include Which samples would you like to include? Usually, you will choose all sample columns, or
#'   you could choose to remove certain samples. Samples excluded here will be removed in this step and from further
#'   analysis downstream of this step. (Default: `NULL` - all sample IDs in `moo@sample_meta` will be used.)
#' @param use_cpm_counts_to_filter If no transformation has been been performed on counts matrix (eg Raw Counts) set to
#'   TRUE. If TRUE counts will be transformed to CPM and filtered based on given criteria. If gene counts matrix has
#'   been transformed (eg log2, CPM, FPKM or some form of Normalization) set to FALSE. If FALSE no further
#'   transformation will be applied and features will be filtered as is. For RNAseq data RAW counts should be
#'   transformed to CPM in order to properly filter.
#' @param minimum_count_value_to_be_considered_nonzero Minimum value in the selected filtering table required for a
#'   sample to be considered nonzero. If `use_cpm_counts_to_filter` is `TRUE`, this threshold is applied to CPM values.
#'   If `use_cpm_counts_to_filter` is `FALSE`, this threshold is applied directly to the selected `count_type` table.
#' @param minimum_number_of_samples_with_nonzero_counts_in_total Minimum number of samples in total that must meet the
#'   `minimum_count_value_to_be_considered_nonzero` threshold for a feature to be kept.
#' @param use_group_based_filtering If TRUE, only keeps features (e.g. genes) that have at least a certain number of
#'   samples passing the threshold in at least one group
#' @param minimum_number_of_samples_with_nonzero_counts_in_a_group Only keeps genes that have at least this number of
#'   samples meeting the threshold in at least one group
#' @param principal_component_on_x_axis The principal component to plot on the x-axis for the PCA plot. Choices include
#'   1, 2, 3, ... (default: 1)
#' @param principal_component_on_y_axis The principal component to plot on the y-axis for the PCA plot. Choices include
#'   1, 2, 3, ... (default: 2)
#' @param legend_position_for_pca legend position for the PCA plot
#' @param point_size_for_pca geom point size for the PCA plot
#' @param add_label_to_pca If `TRUE`, display labels from `label_colname` on PCA points. If `FALSE`, the PCA plot uses
#'   unlabeled points while heatmap labels still use `label_colname`.
#' @param label_font_size label font size for the PCA plot
#' @param label_offset_y_ label offset y for the PCA plot
#' @param label_offset_x_ label offset x for the PCA plot
#' @param samples_to_rename If you do not have a Plot Labels Column in your sample metadata table, you can use this
#'   parameter to rename samples manually for display on the PCA plot. Use "Add item" to add each additional sample for
#'   renaming. Use the following format to describe which old name (in your sample metadata table) you want to rename to
#'   which new name: old_name: new_name
#' @param color_histogram_by_group Set to FALSE to label histogram by Sample Names, or set to TRUE to label histogram by
#'   the column you select in the "Group Column Used to Color Histogram" parameter (below). Default is TRUE.
#' @param set_min_max_for_x_axis_for_histogram whether to set min/max value for histogram x-axis
#' @param minimum_for_x_axis_for_histogram x-axis minimum for histogram plot
#' @param maximum_for_x_axis_for_histogram x-axis maximum for histogram plot
#' @param legend_position_for_histogram legend position for the histogram plot. consider setting to 'none' for a large
#'   number of samples.
#' @param legend_font_size_for_histogram legend font size for the histogram plot.
#'   If `NULL`, the size is scaled automatically.
#' @param number_of_histogram_legend_columns number of columns for the histogram legend
#' @param colors_for_plots Optional colors for PCA/histogram/heatmap plots. If `NULL`, colors are taken from
#'   `moo@analyses$colors[[group_colname]]`.
#'   Colors must either be names in `grDevices::colors()` or valid hex codes.
#'   Unnamed colors are assigned by factor level order when the grouping column is a factor;
#'   otherwise, they follow the order in which groups first appear in the metadata column. If more groups are present
#'   than colors provided,
#'   supplied colors are used first and additional colors are generated from the selected palette for the remaining
#'   groups; random colors are used only if that palette returns fewer colors than the number of groups.
#' @param plot_corr_matrix_heatmap Datasets with a large number of samples may be too large to create a correlation
#'   matrix heatmap. If this function takes longer than 5 minutes to run, Set to `FALSE` and the correlation matrix will
#'   not be be created. Default is `TRUE`.
#' @param interactive_plots set to TRUE to make PCA and Histogram plots interactive with `plotly`, allowing you to hover
#'   your mouse over a point or line to view sample information. The similarity heat map will not display if this toggle
#'   is set to `TRUE`. Default is `FALSE`.
#' @param plots_subdir subdirectory in `figures/` where plots will be saved if `save_plots` is `TRUE`
#'
#' @return `multiOmicDataSet` with filtered counts
#' @export
#'
#' @examples
#' moo <- create_multiOmicDataSet_from_dataframes(
#'   as.data.frame(nidap_sample_metadata),
#'   as.data.frame(nidap_clean_raw_counts),
#'   sample_id_colname = "Sample",
#'   feature_id_colname = "Gene"
#' ) |>
#'   filter_counts(
#'     count_type = "raw"
#'   )
#' head(moo@counts$filt)
#'
#' @family moo methods
#' @family main analysis
filter_counts <- function(
  moo,
  count_type = "clean",
  feature_id_colname = NULL,
  sample_id_colname = NULL,
  group_colname = "Group",
  label_colname = "Label",
  samples_to_include = NULL,
  minimum_count_value_to_be_considered_nonzero = 8,
  minimum_number_of_samples_with_nonzero_counts_in_total = 7,
  minimum_number_of_samples_with_nonzero_counts_in_a_group = 3,
  use_cpm_counts_to_filter = TRUE,
  use_group_based_filtering = FALSE,
  principal_component_on_x_axis = 1,
  principal_component_on_y_axis = 2,
  legend_position_for_pca = "top",
  point_size_for_pca = 5,
  add_label_to_pca = TRUE,
  label_font_size = 3,
  label_offset_y_ = 2,
  label_offset_x_ = 2,
  samples_to_rename = c(""),
  color_histogram_by_group = TRUE,
  set_min_max_for_x_axis_for_histogram = FALSE,
  minimum_for_x_axis_for_histogram = -1,
  maximum_for_x_axis_for_histogram = 1,
  legend_position_for_histogram = "top",
  legend_font_size_for_histogram = NULL,
  number_of_histogram_legend_columns = 6,
  colors_for_plots = NULL,
  plot_corr_matrix_heatmap = TRUE,
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  interactive_plots = FALSE,
  plots_subdir = "filt"
) {
  if (!(count_type %in% names(moo@counts))) {
    stop(glue::glue("count_type {count_type} not in moo@counts"))
  }
  counts_dat <- moo@counts[[count_type]] |> as.data.frame() # currently, this function requires data frames
  sample_metadata <- moo@sample_meta |> as.data.frame()

  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(counts_dat)[1]
  }
  if (is.null(samples_to_include)) {
    samples_to_include <- sample_metadata |> dplyr::pull(sample_id_colname)
  }
  pca_label_colname <- if (isTRUE(add_label_to_pca)) label_colname else NULL
  if (is.null(label_colname)) {
    label_colname <- sample_id_colname
  }
  message(glue::glue("* filtering {count_type} counts"))

  samples_to_include <- samples_to_include |> unlist()

  df <- counts_dat |>
    dplyr::select(
      tidyselect::all_of(feature_id_colname),
      tidyselect::all_of(samples_to_include)
    )

  # filter out low count genes
  df_filt <- remove_low_count_genes(
    counts_dat = df,
    sample_metadata = sample_metadata,
    sample_id_colname = sample_id_colname,
    feature_id_colname = feature_id_colname,
    group_colname = group_colname,
    use_cpm_counts_to_filter = use_cpm_counts_to_filter,
    use_group_based_filtering = use_group_based_filtering,
    minimum_count_value_to_be_considered_nonzero = minimum_count_value_to_be_considered_nonzero,
    minimum_number_of_samples_with_nonzero_counts_in_total = minimum_number_of_samples_with_nonzero_counts_in_total,
    minimum_number_of_samples_with_nonzero_counts_in_a_group = minimum_number_of_samples_with_nonzero_counts_in_a_group
  )
  if (nrow(df_filt) == 0) {
    warning(
      "No features remain after filtering. ",
      "Consider relaxing the filtering thresholds ",
      "(e.g. minimum_count_value_to_be_considered_nonzero, ",
      "minimum_number_of_samples_with_nonzero_counts_in_total, or ",
      "minimum_number_of_samples_with_nonzero_counts_in_a_group). ",
      "Skipping QC plots."
    )
    moo@counts[["filt"]] <- df[0, ]
    return(moo)
  }

  if (isTRUE(print_plots) || isTRUE(save_plots)) {
    default_group_colors <- get_moo_default_colors(
      moo = moo,
      colname = group_colname
    )
    default_label_colors <- get_moo_default_colors(
      moo = moo,
      colname = label_colname
    )
    colors_for_plots <- colors_for_plots %||% default_group_colors

    if (isTRUE(color_histogram_by_group)) {
      colors_for_histogram <- colors_for_plots
    } else {
      colors_for_histogram <- default_label_colors
    }

    pca_plot <- plot_pca(
      df_filt,
      sample_metadata = sample_metadata,
      sample_id_colname = sample_id_colname,
      feature_id_colname = feature_id_colname,
      samples_to_rename = samples_to_rename,
      group_colname = group_colname,
      label_colname = pca_label_colname,
      color_values = colors_for_plots,
      principal_components = c(
        principal_component_on_x_axis,
        principal_component_on_y_axis
      ),
      legend_position = legend_position_for_pca,
      point_size = point_size_for_pca,
      label_font_size = label_font_size,
      label_offset_y_ = label_offset_y_,
      label_offset_x_ = label_offset_x_,
      log_transform = TRUE,
      log_transform_pseudocount = 0.5,
      log_transform_base = "ln",
      print_plots = FALSE,
      save_plots = FALSE
    ) +
      ggplot2::labs(caption = "filtered counts")

    # Histogram data and axis setup
    histogram_x_axis_label <- if (isTRUE(use_cpm_counts_to_filter)) {
      "CPM"
    } else {
      "Count"
    }
    log2_axis_pseudocount <- 0.5
    histogram_threshold_x <- minimum_count_value_to_be_considered_nonzero +
      log2_axis_pseudocount

    # Base histogram plot
    hist_plot <- plot_histogram(
      df_filt,
      sample_metadata,
      sample_id_colname = sample_id_colname,
      feature_id_colname = feature_id_colname,
      group_colname = group_colname,
      label_colname = label_colname,
      color_values = colors_for_histogram,
      color_by_group = color_histogram_by_group,
      set_min_max_for_x_axis = set_min_max_for_x_axis_for_histogram,
      minimum_for_x_axis = minimum_for_x_axis_for_histogram,
      maximum_for_x_axis = maximum_for_x_axis_for_histogram,
      x_axis_label = histogram_x_axis_label,
      legend_position = legend_position_for_histogram,
      legend_font_size = legend_font_size_for_histogram,
      number_of_legend_columns = number_of_histogram_legend_columns,
      interactive_plots = interactive_plots,
      return_ggplot = TRUE,
      use_log2_x_axis = TRUE
    ) +
      ggplot2::labs(caption = "filtered counts")

    # Filtering threshold marker
    hist_plot <- hist_plot +
      ggplot2::geom_vline(
        xintercept = histogram_threshold_x,
        linetype = 2,
        linewidth = 1
      ) +
      ggplot2::annotate(
        "text",
        x = histogram_threshold_x,
        y = Inf,
        label = if (isTRUE(use_cpm_counts_to_filter)) {
          glue::glue("CPM: {minimum_count_value_to_be_considered_nonzero}")
        } else {
          glue::glue("Count: {minimum_count_value_to_be_considered_nonzero}")
        },
        hjust = -0.05,
        vjust = 1.5,
        size = 3
      )

    # Correlation heatmap plot
    if (isTRUE(plot_corr_matrix_heatmap)) {
      corHM <- plot_corr_heatmap(
        df_filt[, samples_to_include],
        sample_metadata = sample_metadata,
        sample_id_colname = sample_id_colname,
        feature_id_colname = feature_id_colname,
        label_colname = label_colname,
        group_colname = group_colname,
        color_values = colors_for_plots
      )
      print_or_save_plot(
        corHM,
        filename = file.path(plots_subdir, "corr_heatmap.png"),
        print_plots = print_plots,
        save_plots = save_plots,
        caption = "filtered counts"
      )
    }

    # Interactive plot conversion
    plot_ext <- "png"
    if (isTRUE(interactive_plots)) {
      pca_plot <- pca_plot |> plotly::ggplotly(tooltip = "text")
      hist_plot <- hist_plot |> plotly::ggplotly(tooltip = "text")
      plot_ext <- "html"
    }

    # Save or print PCA and histogram plots
    if (identical(plot_ext, "png")) {
      print_or_save_plot(
        pca_plot,
        filename = file.path(plots_subdir, glue::glue("pca.{plot_ext}")),
        print_plots = print_plots,
        save_plots = save_plots,
        width = 7,
        height = 7,
        units = "in"
      )
    } else {
      print_or_save_plot(
        pca_plot,
        filename = file.path(plots_subdir, glue::glue("pca.{plot_ext}")),
        print_plots = print_plots,
        save_plots = save_plots
      )
    }
    print_or_save_plot(
      hist_plot,
      filename = file.path(plots_subdir, glue::glue("histogram.{plot_ext}")),
      print_plots = print_plots,
      save_plots = save_plots
    )
  }
  df_final <- df |>
    dplyr::filter(
      !!rlang::sym(feature_id_colname) %in% df_filt[, feature_id_colname]
    )

  moo@counts[["filt"]] <- df_final

  return(moo)
}

#' Remove low-count genes
#'
#' TODO this function also transforms raw counts to CPM, but that should be a separate function before this step, before
#' filter_counts function()
#'
#' @inheritParams filter_counts
#'
#' @return counts matrix with low-count genes removed
#' @keywords internal
#'
remove_low_count_genes <- function(
  counts_dat,
  sample_metadata,
  sample_id_colname = NULL,
  feature_id_colname,
  group_colname,
  use_cpm_counts_to_filter = TRUE,
  use_group_based_filtering = FALSE,
  minimum_count_value_to_be_considered_nonzero = 8,
  minimum_number_of_samples_with_nonzero_counts_in_total = 7,
  minimum_number_of_samples_with_nonzero_counts_in_a_group = 3
) {
  df <- counts_dat

  df <- df[stats::complete.cases(df), ]
  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  sample_colnames <- setdiff(colnames(df), feature_id_colname)
  sample_colnames <- sample_colnames[vapply(
    df[sample_colnames],
    is.numeric,
    logical(1)
  )]

  # USE CPM Transformation
  trans_df <- df
  if (use_cpm_counts_to_filter == TRUE) {
    trans_df[, sample_colnames] <- edgeR::cpm(as.matrix(df[, sample_colnames]))
  }
  passing_counts <- as.matrix(trans_df[, sample_colnames, drop = FALSE]) >=
    minimum_count_value_to_be_considered_nonzero

  if (use_group_based_filtering == TRUE) {
    if (!(sample_id_colname %in% colnames(sample_metadata))) {
      stop(glue::glue(
        "sample_id_colname {sample_id_colname} not in sample_metadata"
      ))
    }
    if (!(group_colname %in% colnames(sample_metadata))) {
      stop(glue::glue("group_colname {group_colname} not in sample_metadata"))
    }
    sample_match <- match(
      sample_colnames,
      as.character(sample_metadata[[sample_id_colname]])
    )
    if (anyNA(sample_match)) {
      missing_samples <- sample_colnames[is.na(sample_match)]
      stop(glue::glue(
        "sample_metadata is missing sample IDs: ",
        "{glue::glue_collapse(missing_samples, sep = ', ')}"
      ))
    }
    sample_groups <- as.character(sample_metadata[[group_colname]][
      sample_match
    ])
    groups <- unique(stats::na.omit(sample_groups))
    if (length(groups) == 0) {
      stop(glue::glue(
        "group_colname {group_colname} has no non-missing values for selected samples"
      ))
    }
    passing_samples_by_group <- vapply(
      groups,
      function(group) {
        group_sample_colnames <- sample_colnames[
          !is.na(sample_groups) & sample_groups == group
        ]
        rowSums(passing_counts[, group_sample_colnames, drop = FALSE])
      },
      numeric(nrow(passing_counts))
    )
    keep_rows <- rowSums(
      passing_samples_by_group >=
        minimum_number_of_samples_with_nonzero_counts_in_a_group
    ) >=
      1
    df_filt <- trans_df[keep_rows, , drop = FALSE]
  } else {
    keep_rows <- (rowSums(passing_counts) >=
      minimum_number_of_samples_with_nonzero_counts_in_total)
    df_filt <- trans_df[keep_rows, , drop = FALSE]
  }
  rownames(df_filt) <- NULL

  message(paste0("Number of features after filtering: ", nrow(df_filt)))
  return(df_filt)
}
