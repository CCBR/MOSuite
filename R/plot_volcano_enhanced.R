#' Enhanced Volcano Plot
#'
#' Uses [Bioconductor's Enhanced Volcano
#' Plot](https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html).
#' An S7 generic with methods for `multiOmicDataSet` and `data.frame`.
#'
#' @param moo_diff multiOmicDataSet or differential expression analysis result data frame.
#'
#' @export
plot_volcano_enhanced <- S7::new_generic(
  "plot_volcano_enhanced",
  "moo_diff",
  function(
    moo_diff,
    feature_id_colname = NULL,
    change_colname = NULL,
    signif_colname = NULL,
    signif_threshold = 0.05,
    change_threshold = 1.0,
    value_to_sort_the_output_dataset = "p-value",
    num_features_to_label = 30,
    label_features = FALSE,
    custom_gene_list = "",
    label_significant_features_only = TRUE,
    label_font_size = 7,
    draw_connectors = FALSE,
    change_sig_name = "p-value",
    change_lfc_name = "log2FC",
    title = "Volcano Plots",
    title_font_size = 24,
    use_custom_lab = FALSE,
    use_default_x_axis_limit = TRUE,
    x_axis_limit = 5,
    use_default_y_axis_limit = TRUE,
    y_axis_limit = 10,
    axis_lab_size = 24,
    axis_tick_lab_size = 16,
    point_size = 2,
    default_label_color = "black",
    custom_label_color = "black",
    color_of_signif_threshold_line = "black",
    color_of_non_significant_features = "grey30",
    color_of_logfold_change_threshold_line = "forestgreen",
    color_of_features_meeting_only_signif_threshold = "royalblue",
    color_for_features_meeting_pvalue_and_foldchange_thresholds = "red2",
    graphics_device = grDevices::png,
    image_width = 3000,
    image_height = 3000,
    dpi = 300,
    use_default_grid_layout = TRUE,
    number_of_rows_in_grid_layout = NULL,
    scale_image_to_grid = FALSE,
    interactive_plots = FALSE,
    print_plots = options::opt("print_plots"),
    save_plots = options::opt("save_plots"),
    plots_subdir = "diff",
    plot_filename = "volcano_enhanced.png"
  ) {
    return(S7::S7_dispatch())
  }
)

#' @rdname plot_volcano_enhanced
S7::method(plot_volcano_enhanced, multiOmicDataSet) <- function(
  moo_diff,
  feature_id_colname = NULL,
  change_colname = NULL,
  signif_colname = NULL,
  signif_threshold = 0.05,
  change_threshold = 1.0,
  value_to_sort_the_output_dataset = "p-value",
  num_features_to_label = 30,
  label_features = FALSE,
  custom_gene_list = "",
  label_significant_features_only = TRUE,
  label_font_size = 7,
  draw_connectors = FALSE,
  change_sig_name = "p-value",
  change_lfc_name = "log2FC",
  title = "Volcano Plots",
  title_font_size = 24,
  use_custom_lab = FALSE,
  use_default_x_axis_limit = TRUE,
  x_axis_limit = 5,
  use_default_y_axis_limit = TRUE,
  y_axis_limit = 10,
  axis_lab_size = 24,
  axis_tick_lab_size = 16,
  point_size = 2,
  default_label_color = "black",
  custom_label_color = "black",
  color_of_signif_threshold_line = "black",
  color_of_non_significant_features = "grey30",
  color_of_logfold_change_threshold_line = "forestgreen",
  color_of_features_meeting_only_signif_threshold = "royalblue",
  color_for_features_meeting_pvalue_and_foldchange_thresholds = "red2",
  graphics_device = grDevices::png,
  image_width = 3000,
  image_height = 3000,
  dpi = 300,
  use_default_grid_layout = TRUE,
  number_of_rows_in_grid_layout = NULL,
  scale_image_to_grid = FALSE,
  interactive_plots = FALSE,
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = "diff",
  plot_filename = "volcano_enhanced.png"
) {
  return(
    join_dfs_wide(moo_diff@analyses$diff) |>
      plot_volcano_enhanced(
        feature_id_colname,
        change_colname,
        signif_colname,
        signif_threshold,
        change_threshold,
        value_to_sort_the_output_dataset,
        num_features_to_label,
        label_features = label_features,
        custom_gene_list = custom_gene_list,
        label_significant_features_only = label_significant_features_only,
        label_font_size = label_font_size,
        draw_connectors = draw_connectors,
        change_sig_name = change_sig_name,
        change_lfc_name = change_lfc_name,
        title = title,
        title_font_size = title_font_size,
        use_custom_lab = use_custom_lab,
        use_default_x_axis_limit = use_default_x_axis_limit,
        x_axis_limit = x_axis_limit,
        use_default_y_axis_limit = use_default_y_axis_limit,
        y_axis_limit = y_axis_limit,
        axis_lab_size = axis_lab_size,
        axis_tick_lab_size = axis_tick_lab_size,
        point_size = point_size,
        default_label_color = default_label_color,
        custom_label_color = custom_label_color,
        color_of_signif_threshold_line = color_of_signif_threshold_line,
        color_of_non_significant_features = color_of_non_significant_features,
        color_of_logfold_change_threshold_line = color_of_logfold_change_threshold_line,
        color_of_features_meeting_only_signif_threshold = color_of_features_meeting_only_signif_threshold,
        color_for_features_meeting_pvalue_and_foldchange_thresholds = color_for_features_meeting_pvalue_and_foldchange_thresholds,
        graphics_device = graphics_device,
        image_width = image_width,
        image_height = image_height,
        dpi = dpi,
        use_default_grid_layout = use_default_grid_layout,
        number_of_rows_in_grid_layout = number_of_rows_in_grid_layout,
        scale_image_to_grid = scale_image_to_grid,
        interactive_plots = interactive_plots,
        print_plots = print_plots,
        save_plots = save_plots,
        plots_subdir = plots_subdir,
        plot_filename = plot_filename
      )
  )
}

#' @inheritParams option_params
#' @inheritParams filter_counts
#'
#' @param moo_diff Differential expression analysis result from one or more contrasts. This must be a dataframe.
#' @param change_colname Character vector of logFC column names, one per
#'   contrast (e.g. `c("B-A_logFC", "C-A_logFC")`). Defaults to `NULL`, which
#'   auto-detects all columns ending in `_logFC`.
#' @param signif_colname Character vector of significance column names, one per
#'   contrast (e.g. `c("B-A_adjpval", "C-A_adjpval")`). Defaults to `NULL`,
#'   which auto-detects corresponding columns by checking for `_adjpval` first,
#'   then `_pval`, for each contrast in `change_colname`.
#' @param signif_threshold Numeric significance threshold (p-value or adjusted p-value cutoff). Default: 0.05
#' @param change_threshold Numeric value specifying the fold change cutoff for significance (i.e. filters on
#'   `change_colname`)
#' @param value_to_sort_the_output_dataset How to sort the output dataset. Options are "fold-change", "p-value", or
#'   "t-statistic".
#' @param num_features_to_label Number of top features/genes to label in the volcano plot. Default is 30.
#' @param label_features If `TRUE`, only the features specified in `custom_gene_list` will be used for labeling in the
#'   volcano plot, ignoring the top features.
#' @param custom_gene_list comma-separated string of feature names or IDs to include in the volcano plot.
#' @param label_significant_features_only If `TRUE`, automatic labels are selected only from features that pass both
#'   the significance and fold-change thresholds.
#' @param label_font_size Size of the labels in the volcano plot.
#' @param draw_connectors If `TRUE`, draw connector lines from labels to their points and spread labels to reduce
#'   overlap.
#' @param change_sig_name Name for the significance column in the plot. Default is "p-value".
#' @param change_lfc_name Name for the fold change column in the plot. Default is "log2FC".
#' @param title Title of the plot. Default is "Volcano Plots".
#' @param title_font_size Size of the plot title.
#' @param use_custom_lab If TRUE, uses custom labels for the plot (set by `change_sig_name` and `change_lfc_name`)
#' @param use_default_x_axis_limit Set to TRUE to use the default x-axis limit.
#' @param x_axis_limit Custom x-axis limit. A single value is treated symmetrically, and a two-value vector is treated
#'   as lower and upper limits.
#' @param use_default_y_axis_limit Set to TRUE to use the default y-axis limit.
#' @param y_axis_limit Custom y-axis limit.
#' @param axis_lab_size Size of the axis labels.
#' @param axis_tick_lab_size Size of the axis tick labels.
#' @param point_size Size of the points in the plot.
#' @param default_label_color Set the color for the text used to add feature labels to points.
#' @param custom_label_color Set the color for labels from `custom_gene_list`.
#' @param color_of_signif_threshold_line Color of the significance threshold line.
#' @param color_of_non_significant_features Color of the non-significant features.
#' @param color_of_logfold_change_threshold_line Color of the features that meet only the log fold change threshold.
#' @param color_of_features_meeting_only_signif_threshold Color of the features that meet only the significance
#'   threshold.
#' @param color_for_features_meeting_pvalue_and_foldchange_thresholds Color of the features that meet both the p-value
#'   and fold change thresholds.
#' @param graphics_device passed to `ggsave(device)`. Default: `grDevices::png`
#' @param image_width output image width in pixels - only used if save_plots is TRUE
#' @param image_height output image height in pixels - only used if save_plots is TRUE
#' @param dpi dots-per-inch of the output image (see `ggsave()`) - only used if save_plots is TRUE
#' @param use_default_grid_layout Retained for compatibility. Grid layout is handled by `plot_volcano_summary()`.
#' @param number_of_rows_in_grid_layout Retained for compatibility. Grid layout is handled by `plot_volcano_summary()`.
#' @param plot_filename plot output filename - only used if save_plots is TRUE. When multiple comparisons are saved
#'   separately, the comparison name is appended before the file extension.
#' @param scale_image_to_grid Retained for compatibility. Grid layout is handled by `plot_volcano_summary()`.
#'
#' @keywords plotters volcano
#'
#' @examples
#' plot_volcano_enhanced(nidap_deg_analysis, print_plots = TRUE)
#'
#' @rdname plot_volcano_enhanced
S7::method(plot_volcano_enhanced, S7::class_data.frame) <- function(
  moo_diff,
  feature_id_colname = NULL,
  change_colname = NULL,
  signif_colname = NULL,
  signif_threshold = 0.05,
  change_threshold = 1.0,
  value_to_sort_the_output_dataset = "p-value",
  num_features_to_label = 30,
  label_features = FALSE,
  custom_gene_list = "",
  label_significant_features_only = TRUE,
  label_font_size = 7,
  draw_connectors = FALSE,
  change_sig_name = "p-value",
  change_lfc_name = "log2FC",
  title = "Volcano Plots",
  title_font_size = 24,
  use_custom_lab = FALSE,
  use_default_x_axis_limit = TRUE,
  x_axis_limit = 5,
  use_default_y_axis_limit = TRUE,
  y_axis_limit = 10,
  axis_lab_size = 24,
  axis_tick_lab_size = 16,
  point_size = 2,
  default_label_color = "black",
  custom_label_color = "black",
  color_of_signif_threshold_line = "black",
  color_of_non_significant_features = "grey30",
  color_of_logfold_change_threshold_line = "forestgreen",
  color_of_features_meeting_only_signif_threshold = "royalblue",
  color_for_features_meeting_pvalue_and_foldchange_thresholds = "red2",
  graphics_device = grDevices::png,
  image_width = 3000,
  image_height = 3000,
  dpi = 300,
  use_default_grid_layout = TRUE,
  number_of_rows_in_grid_layout = NULL,
  scale_image_to_grid = FALSE,
  interactive_plots = FALSE,
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = "diff",
  plot_filename = "volcano_enhanced.png"
) {
  abort_packages_not_installed("EnhancedVolcano")
  ### PH
  # Input - DEG table from Limma DEG template
  # Output - Volcano plot + interactive Volcano Plot
  # Purpose Create detailed Volcano for each contrast individually

  diff_dat <- as.data.frame(moo_diff)

  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(diff_dat)[1]
  }
  label_col <- feature_id_colname

  # Resolve/auto-detect change_colname and signif_colname
  resolved <- resolve_volcano_colnames(diff_dat, change_colname, signif_colname)
  change_colname <- resolved$change_colname
  signif_colname <- resolved$signif_colname

  rank <- list()
  plots_list <- list()

  # user can select multiple comparisons to create volcano plots
  for (i in seq_along(change_colname)) {
    ### PH: START Build table for Volcano plot

    lfccol <- change_colname[i]
    sigcol <- signif_colname[i]
    columns_of_interest <- c(label_col, change_colname[i], signif_colname[i])
    df <- diff_dat |>
      dplyr::select(tidyselect::one_of(columns_of_interest)) |>
      dplyr::mutate(
        !!rlang::sym(lfccol) := tidyr::replace_na(!!rlang::sym(lfccol), 0)
      ) |>
      dplyr::mutate(
        !!rlang::sym(sigcol) := tidyr::replace_na(!!rlang::sym(sigcol), 1)
      )
    # mutate(.data[[lfc.col[i]]] = replace_na(.data[[lfc.col[i]]], 0)) |>
    # mutate(.data[[sig.col[i]]] = replace_na(.data[[sig.col[i]]], 1))
    if (use_custom_lab == TRUE) {
      lfc_name <- if (nchar(change_lfc_name) == 0) {
        change_colname[i]
      } else {
        change_lfc_name
      }
      sig_name <- if (nchar(change_sig_name) == 0) {
        signif_colname[i]
      } else {
        change_sig_name
      }
      colnames(df) <- c(label_col, lfc_name, sig_name)
    } else {
      lfc_name <- change_colname[i]
      sig_name <- signif_colname[i]
    }

    ### PH: START Creating rank based on pvalue and fold change
    # This is unique to this template and could be useful as a generic tool to create rankes for GSEA. Recommend
    # extracting this function
    group <- gsub("_pval|p_val_", "", sig_name)
    rank[[i]] <- -log10(df[[sig_name]]) * sign(df[[lfc_name]])
    names(rank)[i] <- paste0("C_", group, "_rank")
    ### PH: End Creating rank based on pvalue and fold change

    message(paste0("Genes in initial dataset: ", nrow(df), "\n"))

    # Select top genes by logFC or Significance
    contrast_label <- gsub("_logFC$", "", change_colname[i])
    tstat_colname <- paste0(contrast_label, "_tstat")
    if (value_to_sort_the_output_dataset == "fold-change") {
      df <- df |> dplyr::arrange(dplyr::desc(.data[[lfc_name]]))
    } else if (value_to_sort_the_output_dataset == "p-value") {
      df <- df |> dplyr::arrange(.data[[sig_name]])
    } else if (value_to_sort_the_output_dataset == "t-statistic") {
      if (tstat_colname %in% colnames(diff_dat)) {
        df <- df |>
          dplyr::mutate(.mosuite_sort_tstat = diff_dat[[tstat_colname]]) |>
          dplyr::arrange(dplyr::desc(abs(.data$.mosuite_sort_tstat))) |>
          # Previous tidyselect form: dplyr::select(-.data$.mosuite_sort_tstat)
          dplyr::select(-tidyselect::all_of(".mosuite_sort_tstat"))
      } else {
        warning(glue::glue(
          "Could not find t-statistic column '{tstat_colname}'. Labels were not sorted by t-statistic."
        ))
      }
    }

    if (label_significant_features_only) {
      df_sub <- df[
        df[[sig_name]] <= signif_threshold &
          abs(df[[lfc_name]]) >= change_threshold,
      ]
    } else {
      df_sub <- df
    }

    genes_to_label <- as.character(df_sub[1:num_features_to_label, label_col])
    split_values <- unlist(strsplit(gsub(",", " ", custom_gene_list), " "))
    custom_labels <- split_values[split_values != ""]

    filter <- custom_labels %in% df[, label_col]
    missing_labels <- custom_labels[!filter]
    custom_labels <- custom_labels[filter]

    if (length(missing_labels) > 0) {
      message(glue::glue(
        ("Could not find missing labels:\t{paste(missing_labels, collapse = ', ')}")
      ))
    }

    if (label_features) {
      genes_to_label <- custom_labels
    } else {
      genes_to_label <- unique(append(genes_to_label, custom_labels))
    }

    labels_in_plot_order <- df[[label_col]][df[[label_col]] %in% genes_to_label]
    label_colors <- rep(default_label_color, length(labels_in_plot_order))
    label_colors[labels_in_plot_order %in% custom_labels] <- custom_label_color
    if (length(label_colors) == 0) {
      label_colors <- default_label_color
    }

    significant <- vector(length = nrow(df))
    significant[] <- "Not significant"
    significant[which(abs(df[, 2]) > change_threshold)] <- "Fold change only"
    significant[which(df[, 3] < signif_threshold)] <- "Significant only"
    significant[which(
      abs(df[, 2]) > change_threshold & df[, 3] < signif_threshold
    )] <- "Significant and fold change"

    ### PH: END Build table for Volcano plot

    ### PH: START Create Volcano plot

    ### PH: Set Axis limits - Unique feature to this plot that should be included with any Volcano plot function
    ##############################

    ## Y-axis range change:
    # fix pvalue == 0
    shapeCustom <- rep(19, nrow(df))
    maxy <- max(-log10(df[[sig_name]]), na.rm = TRUE)
    if (!isTRUE(use_default_y_axis_limit)) {
      maxy <- y_axis_limit
    }

    message(paste0("Max y: ", maxy, "\n"))
    if (maxy == Inf) {
      # Sometimes, pvalues == 0
      keep <- df[[sig_name]] > 0
      df[[sig_name]][!keep] <- min(df[[sig_name]][keep])
      shapeCustom[!keep] <- 17

      maxy <- -log10(min(df[[sig_name]][keep]))
      message("Some p-values equal zero. Adjusting y-limits.\n")
      message(paste0("Max y adjusted: ", maxy, "\n"))
    }

    # By default, nothing will be greater than maxy. User can set this value lower
    keep <- -log10(df[[sig_name]]) <= maxy
    df[[sig_name]][!keep] <- maxy
    shapeCustom[!keep] <- 17

    names(shapeCustom) <- rep("Exact", length(shapeCustom))
    names(shapeCustom)[shapeCustom == 17] <- "Adjusted"

    # Remove if nothin' doin'
    if (all(shapeCustom == 19)) {
      shapeCustom <- NULL
    }
    maxy <- ceiling(maxy)

    ## X-axis custom range change:
    if (isTRUE(use_default_x_axis_limit)) {
      xlim <- c(
        floor(min(df[, lfc_name])),
        ceiling(max(df[, lfc_name]))
      )
    } else if (is.character(x_axis_limit) && grepl(",", x_axis_limit)) {
      split_values <- strsplit(x_axis_limit, ",")[[1]]
      xlim <- as.numeric(trimws(split_values[1:2]))
    } else if (length(x_axis_limit) == 1) {
      xlim <- c(
        -1 * as.numeric(x_axis_limit),
        as.numeric(x_axis_limit)
      )
    } else {
      xlim <- as.numeric(x_axis_limit[1:2])
    }

    ### Create axis labels
    ##############################

    if (grepl("log", change_colname[i])) {
      xlab <- bquote(~ Log[2] ~ "fold change")
    } else {
      xlab <- "Fold change"
    }
    if (grepl("adj", signif_colname[i])) {
      ylab <- bquote(~ -Log[10] ~ "FDR")
      signif_legend_label <- "FDR"
      signif_and_fc_legend_label <- "FDR and log2 FC"
    } else {
      ylab <- bquote(~ -Log[10] ~ "p-value")
      signif_legend_label <- "p-value"
      signif_and_fc_legend_label <- "p-value and log2 FC"
    }
    legend_labels <- c(
      "NS",
      expression(Log[2] ~ FC),
      signif_legend_label,
      signif_and_fc_legend_label
    )
    if (use_custom_lab) {
      if (lfc_name != change_colname[i]) {
        xlab <- gsub("_", " ", lfc_name)
      }
      if (sig_name != signif_colname[i]) {
        ylab <- gsub("_", " ", sig_name)
      }
    }

    title_is_default <- identical(title, "Volcano Plots")
    title_is_per_plot <- length(title) == length(change_colname)
    plot_title <- if (title_is_default) {
      contrast_label
    } else if (title_is_per_plot) {
      title[i]
    } else {
      title
    }
    plot_subtitle <- if (title_is_default || title_is_per_plot) NULL else group

    volcano_plot <- EnhancedVolcano::EnhancedVolcano(
      df,
      x = lfc_name,
      y = sig_name,
      lab = df[, label_col],
      selectLab = genes_to_label,
      title = plot_title,
      subtitle = plot_subtitle,
      xlab = xlab,
      ylab = ylab,
      xlim = xlim,
      ylim = c(0, maxy),
      pCutoff = signif_threshold,
      FCcutoff = change_threshold,
      axisLabSize = axis_lab_size,
      titleLabSize = title_font_size,
      legendLabels = legend_labels,
      labSize = label_font_size,
      drawConnectors = draw_connectors,
      labCol = label_colors,
      pointSize = point_size,
      col = c(
        color_of_non_significant_features,
        color_of_logfold_change_threshold_line,
        color_of_features_meeting_only_signif_threshold,
        color_for_features_meeting_pvalue_and_foldchange_thresholds
      ),
      cutoffLineCol = color_of_signif_threshold_line,
      shapeCustom = shapeCustom
    ) +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = axis_tick_lab_size)
      )

    ## Creating plot that can be converted to plotly interactive plot (no labels):
    ## PH: make this feature an option not default
    if (isTRUE(interactive_plots)) {
      p_empty <- EnhancedVolcano::EnhancedVolcano(
        df,
        x = lfc_name,
        y = sig_name,
        lab = rep("", nrow(df)),
        # Setting labels to empty strings
        selectLab = NULL,
        title = plot_title,
        subtitle = plot_subtitle,
        xlab = xlab,
        ylab = ylab,
        xlim = xlim,
        ylim = c(0, maxy),
        pCutoff = signif_threshold,
        FCcutoff = change_threshold,
        axisLabSize = axis_lab_size,
        titleLabSize = title_font_size,
        legendLabels = legend_labels,
        labSize = label_font_size,
        drawConnectors = draw_connectors,
        labCol = default_label_color,
        pointSize = point_size,
        col = c(
          color_of_non_significant_features,
          color_of_logfold_change_threshold_line,
          color_of_features_meeting_only_signif_threshold,
          color_for_features_meeting_pvalue_and_foldchange_thresholds
        ),
        cutoffLineCol = color_of_signif_threshold_line,
        shapeCustom = shapeCustom
      ) +
        ggplot2::theme(
          axis.text = ggplot2::element_text(size = axis_tick_lab_size)
        )

      # Extract the data used for plotting
      plot_data <- ggplot2::ggplot_build(p_empty)$data[[1]]

      pxx <- p_empty +
        ggplot2::xlab("Fold Change") + # Simplify x-axis label
        ggplot2::ylab("Significance") + # Simplify y-axis label
        ggplot2::theme_minimal() +
        ggplot2::geom_point(
          ggplot2::aes(
            text = paste(
              "Gene:",
              df[[label_col]],
              "<br>change_threshold:",
              df[[lfc_name]],
              "<br>P-value:",
              df[[sig_name]]
            ),
            colour = as.character(plot_data$colour),
            fill = as.character(plot_data$colour) # Set fill to the same as colour
          ),
          shape = 21,
          # Shape that supports both colour and fill
          size = 2,
          # Size of the points
          stroke = 0.1 # Stroke width
        ) +
        ggplot2::scale_fill_identity()

      # Add interactive hover labels for the gene names
      volcano_plot <- plotly::ggplotly(pxx, tooltip = c("text"))
    }
    plots_list[[i]] <- volcano_plot
  }
  suffix_plot_filename <- function(filename, suffix) {
    filename_dir <- dirname(filename)
    filename_base <- basename(filename)
    extension <- tools::file_ext(filename_base)
    filename_stem <- if (nzchar(extension)) {
      sub(paste0("\\.", extension, "$"), "", filename_base)
    } else {
      filename_base
    }
    safe_suffix <- gsub("[^[:alnum:]_.-]+", "_", suffix)
    safe_suffix <- gsub("^_+|_+$", "", safe_suffix)
    if (!nzchar(safe_suffix)) {
      safe_suffix <- "plot"
    }

    suffixed_filename <- paste0(
      filename_stem,
      "_",
      safe_suffix,
      if (nzchar(extension)) paste0(".", extension) else ""
    )
    if (filename_dir %in% c(".", "")) {
      suffixed_filename
    } else {
      file.path(filename_dir, suffixed_filename)
    }
  }

  nplots <- length(plots_list)
  plot_output_filenames <- rep(plot_filename, nplots)
  if (nplots > 1) {
    comparison_labels <- gsub("_logFC$", "", change_colname[seq_len(nplots)])
    plot_output_filenames <- vapply(
      comparison_labels,
      function(comparison_label) {
        suffix_plot_filename(plot_filename, comparison_label)
      },
      character(1)
    )
  }

  for (i in seq_along(plots_list)) {
    print_or_save_plot(
      plots_list[[i]],
      filename = file.path(plots_subdir, plot_output_filenames[[i]]),
      print_plots = print_plots,
      save_plots = save_plots,
      units = "px",
      width = image_width,
      height = image_height,
      dpi = dpi,
      device = graphics_device
    )
  }

  df_final <- cbind(diff_dat, do.call(cbind, rank))
  attr(df_final, "plots") <- plots_list
  return(df_final)
}
