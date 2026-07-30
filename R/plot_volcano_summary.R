#' Volcano Plot - Summary
#'
#' Produces one volcano plot for each tested contrast in the input DEG table.
#' It can be sorted by either fold change, t-statistic, or p-value. The returned dataset includes one row for each
#' significant gene in each contrast, and contains columns from the DEG analysis of that contrast as well as columns
#' useful to the Venn diagram template downstream.
#' An S7 generic with methods for `multiOmicDataSet` and `data.frame`.
#'
#' @param moo_diff multiOmicDataSet or differential expression analysis result data frame.
#'
#' @export
plot_volcano_summary <- S7::new_generic(
  "plot_volcano_summary",
  "moo_diff",
  function(
    moo_diff,
    feature_id_colname = NULL,
    change_colname = NULL,
    signif_colname = NULL,
    signif_threshold = 0.05,
    change_threshold = 1,
    value_to_sort_the_output_dataset = "t-statistic",
    num_features_to_label = 30,
    add_features = FALSE,
    label_features = FALSE,
    custom_gene_list = "",
    label_significant_features_only = TRUE,
    default_label_color = "black",
    custom_label_color = "black",
    label_font_size = 7,
    draw_connectors = FALSE,
    change_sig_name = "p-value",
    change_lfc_name = "log2FC",
    title = "Volcano Plots",
    title_font_size = 24,
    use_custom_lab = FALSE,
    color_of_signif_threshold_line = "black",
    color_of_non_significant_features = "grey30",
    color_of_logfold_change_threshold_line = "forestgreen",
    color_of_features_meeting_only_signif_threshold = "royalblue",
    color_for_features_meeting_pvalue_and_foldchange_thresholds = "red2",
    flip_vplot = FALSE,
    use_default_x_axis_limit = TRUE,
    x_axis_limit = 5,
    use_default_y_axis_limit = TRUE,
    y_axis_limit = 10,
    point_size = 2,
    axis_lab_size = 24,
    axis_tick_lab_size = 16,
    add_deg_columns = c("FC", "logFC", "tstat", "pval", "adjpval"),
    graphics_device = grDevices::png,
    image_width = 15,
    image_height = 15,
    dpi = 300,
    use_default_grid_layout = TRUE,
    number_of_rows_in_grid_layout = 1,
    plot_filename = "volcano_summary.png",
    print_plots = options::opt("print_plots"),
    save_plots = options::opt("save_plots"),
    plots_subdir = "diff"
  ) {
    return(S7::S7_dispatch())
  }
)

#' @rdname plot_volcano_summary
S7::method(plot_volcano_summary, multiOmicDataSet) <- function(
  moo_diff,
  feature_id_colname = NULL,
  change_colname = NULL,
  signif_colname = NULL,
  signif_threshold = 0.05,
  change_threshold = 1,
  value_to_sort_the_output_dataset = "t-statistic",
  num_features_to_label = 30,
  add_features = FALSE,
  label_features = FALSE,
  custom_gene_list = "",
  label_significant_features_only = TRUE,
  default_label_color = "black",
  custom_label_color = "black",
  label_font_size = 7,
  draw_connectors = FALSE,
  change_sig_name = "p-value",
  change_lfc_name = "log2FC",
  title = "Volcano Plots",
  title_font_size = 24,
  use_custom_lab = FALSE,
  color_of_signif_threshold_line = "black",
  color_of_non_significant_features = "grey30",
  color_of_logfold_change_threshold_line = "forestgreen",
  color_of_features_meeting_only_signif_threshold = "royalblue",
  color_for_features_meeting_pvalue_and_foldchange_thresholds = "red2",
  flip_vplot = FALSE,
  use_default_x_axis_limit = TRUE,
  x_axis_limit = 5,
  use_default_y_axis_limit = TRUE,
  y_axis_limit = 10,
  point_size = 2,
  axis_lab_size = 24,
  axis_tick_lab_size = 16,
  add_deg_columns = c("FC", "logFC", "tstat", "pval", "adjpval"),
  graphics_device = grDevices::png,
  image_width = 15,
  image_height = 15,
  dpi = 300,
  use_default_grid_layout = TRUE,
  number_of_rows_in_grid_layout = 1,
  plot_filename = "volcano_summary.png",
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = "diff"
) {
  return(
    moo_diff@analyses$diff |>
      join_dfs_wide() |>
      plot_volcano_summary(
        feature_id_colname,
        change_colname,
        signif_colname,
        signif_threshold,
        change_threshold,
        value_to_sort_the_output_dataset,
        num_features_to_label,
        add_features,
        label_features,
        custom_gene_list,
        label_significant_features_only,
        default_label_color,
        custom_label_color,
        label_font_size,
        draw_connectors,
        change_sig_name,
        change_lfc_name,
        title,
        title_font_size,
        use_custom_lab,
        color_of_signif_threshold_line,
        color_of_non_significant_features,
        color_of_logfold_change_threshold_line,
        color_of_features_meeting_only_signif_threshold,
        color_for_features_meeting_pvalue_and_foldchange_thresholds,
        flip_vplot,
        use_default_x_axis_limit,
        x_axis_limit,
        use_default_y_axis_limit,
        y_axis_limit,
        point_size,
        axis_lab_size,
        axis_tick_lab_size,
        add_deg_columns,
        graphics_device,
        image_width,
        image_height,
        dpi,
        use_default_grid_layout,
        number_of_rows_in_grid_layout,
        plot_filename,
        print_plots,
        save_plots,
        plots_subdir
      )
  )
}

#' @inheritParams option_params
#' @inheritParams plot_volcano_enhanced
#' @inheritParams filter_counts
#'
#' @param change_colname Character vector of full logFC column names, one per
#'   contrast (e.g. `c("B-A_logFC", "C-A_logFC")`). Defaults to `NULL`, which
#'   auto-detects all columns ending in `_logFC`.
#' @param signif_colname Character vector of full significance column names, one
#'   per contrast (e.g. `c("B-A_adjpval", "C-A_adjpval")`). Defaults to `NULL`,
#'   which auto-detects corresponding columns by checking for `_adjpval` first,
#'   then `_pval`, for each contrast in `change_colname`.
#' @param signif_threshold Numeric significance threshold (p-value or adjusted p-value cutoff). Default: 0.05
#' @param add_features Add custom_gene_list To Labels. Set TRUE when you want to label a specific set of features
#'   (features) in the "custom_gene_list" parameter" IN ADDITION to the number of features you set in the "Number of
#'   Features to Label" parameter.
#' @param label_features Select TRUE when you want to label ONLY a specific list of features(features) given in the
#'   "custom_gene_list" parameter.
#' @param custom_gene_list Provide a list of features (comma separated) to be labeled on the volcano plot. You must
#'   toggle one of the following ON to see these labels: "Add features" or "Label Only My Feature List".
#' @param label_significant_features_only If `TRUE`, automatic labels are selected only from features that pass both
#'   the significance and fold-change thresholds.
#' @param default_label_color Set the color for the text used to add feature (gene) name labels to points.
#' @param custom_label_color Set the color for the specific list of features (features) provided in the "Feature List"
#'   parameter.
#' @param label_font_size Set the font size of the labels. Default: 7
#' @param draw_connectors If `TRUE`, draw connector lines from labels to their points and spread labels to reduce
#'   overlap.
#' @param change_sig_name Name for the significance column in the plot. Default is "p-value".
#' @param change_lfc_name Name for the fold change column in the plot. Default is "log2FC".
#' @param title Title of the plot. Default is "Volcano Plots".
#' @param title_font_size Size of the plot title. Default: 24
#' @param use_custom_lab If TRUE, uses custom labels for the plot axes, set by `change_sig_name` and
#'   `change_lfc_name`.
#' @param color_of_signif_threshold_line Color of the significance threshold line. Default: "black"
#' @param color_of_non_significant_features Color of the non-significant features. Default: "grey30"
#' @param color_of_logfold_change_threshold_line Color of the log fold change threshold line. Default: "forestgreen"
#' @param color_of_features_meeting_only_signif_threshold Color of the features that meet only the significance
#'   threshold. Default: "royalblue"
#' @param color_for_features_meeting_pvalue_and_foldchange_thresholds Color of the features that meet both the p-value
#'   and fold change thresholds. Default: "red2"
#' @param flip_vplot Set to TRUE to flip the fold change values so that the volcano plot looks like a comparison was
#'   B-A. Default: FALSE
#' @param use_default_x_axis_limit Set to TRUE to use the default x-axis limit. Default: TRUE
#' @param x_axis_limit Custom x-axis limit. Default: c(-5, 5)
#' @param use_default_y_axis_limit Set to TRUE to use the default y-axis limit. Default: TRUE
#' @param y_axis_limit Custom y-axis limit. Default: c(0, 10)
#' @param point_size Size of the points in the plot. Default: 1
#' @param axis_lab_size Size of the axis labels. Default: 24
#' @param axis_tick_lab_size Size of the axis tick labels. Default: 16
#' @param add_deg_columns Add additional columns from the DEG analysis to the
#'   output dataset. Default: `"FC", "logFC", "tstat", "pval", "adjpval"`
#' @param use_default_grid_layout Set to TRUE to use the default grid layout. Default: TRUE
#' @param number_of_rows_in_grid_layout Number of rows in the grid layout. Default: 1
#' @param graphics_device passed to `ggsave(device)`. Default: `grDevices::png`
#' @param plot_filename Filename for the output plot. Default: "volcano_plot.png"
#'
#' @keywords plotters volcano
#'
#' @examples
#' plot_volcano_summary(nidap_deg_analysis, print_plots = TRUE)
#'
#' @rdname plot_volcano_summary
#'
S7::method(plot_volcano_summary, S7::class_data.frame) <- function(
  moo_diff,
  feature_id_colname = NULL,
  change_colname = NULL,
  signif_colname = NULL,
  signif_threshold = 0.05,
  change_threshold = 1,
  value_to_sort_the_output_dataset = "t-statistic",
  num_features_to_label = 30,
  add_features = FALSE,
  label_features = FALSE,
  custom_gene_list = "",
  label_significant_features_only = TRUE,
  default_label_color = "black",
  custom_label_color = "black",
  label_font_size = 7,
  draw_connectors = FALSE,
  change_sig_name = "p-value",
  change_lfc_name = "log2FC",
  title = "Volcano Plots",
  title_font_size = 24,
  use_custom_lab = FALSE,
  color_of_signif_threshold_line = "black",
  color_of_non_significant_features = "grey30",
  color_of_logfold_change_threshold_line = "forestgreen",
  color_of_features_meeting_only_signif_threshold = "royalblue",
  color_for_features_meeting_pvalue_and_foldchange_thresholds = "red2",
  flip_vplot = FALSE,
  use_default_x_axis_limit = TRUE,
  x_axis_limit = 5,
  use_default_y_axis_limit = TRUE,
  y_axis_limit = 10,
  point_size = 2,
  axis_lab_size = 24,
  axis_tick_lab_size = 16,
  add_deg_columns = c("FC", "logFC", "tstat", "pval", "adjpval"),
  graphics_device = grDevices::png,
  image_width = 15,
  image_height = 15,
  dpi = 300,
  use_default_grid_layout = TRUE,
  number_of_rows_in_grid_layout = 1,
  plot_filename = "volcano_summary.png",
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = "diff"
) {
  abort_packages_not_installed("EnhancedVolcano")
  diff_dat <- as.data.frame(moo_diff)

  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##
  ### PH
  # Input - DEG table from Limma DEG template
  # Output - Venn Diagrams for selected Comparisons +
  #     Simplified DEG table for selected Comparisons (Only used for Venn Diagram)
  # Purpose - Create Multiple Venn Diagrams
  # Can we use Visualizations from Advanced Volcano function used in the stand alone Volcano plot?

  # Basic Parameters:
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(diff_dat)[1]
  }

  #  Identify all contrasts in DEG output table
  resolved <- resolve_volcano_colnames(diff_dat, change_colname, signif_colname)
  change_colname <- resolved$change_colname
  signif_colname <- resolved$signif_colname
  contrasts <- unique(sub("_logFC$", "", change_colname))

  df_outs <- list()
  plot_change_colnames <- character(0)
  plot_signif_colnames <- character(0)

  #  Create Volcano for each DEG comparison
  for (i in seq_along(contrasts)) {
    contrast <- contrasts[i]
    ### PH: START Build table for Volcano plot
    message(paste0("Preparing table for contrast: ", contrast))
    lfccol <- change_colname[i]
    pvalcol <- signif_colname[i]
    tstatcol <- paste0(contrast, "_", "tstat")

    message(paste0("Fold change column: ", lfccol))
    message(paste0("Significance column: ", pvalcol))
    message(paste0(
      "Total number of features included in volcano plot: ",
      nrow(diff_dat)
    ))

    if (value_to_sort_the_output_dataset == "fold-change") {
      diff_dat <- diff_dat |>
        dplyr::arrange(dplyr::desc(abs(diff_dat[, lfccol])))
    } else if (value_to_sort_the_output_dataset == "p-value") {
      diff_dat <- diff_dat |> dplyr::arrange(diff_dat[, pvalcol])
    } else if (value_to_sort_the_output_dataset == "t-statistic") {
      diff_dat <- diff_dat |>
        dplyr::arrange(dplyr::desc(abs(diff_dat[, tstatcol])))
    }
    feature_ids <- diff_dat[[feature_id_colname]]

    ## optional Parameter: IF DEG was set up A-B User can Flip FC values so that Volcano plot looks like comparison was
    ## B-A
    ## flip contrast section
    indc <- which(colnames(diff_dat) == lfccol) # get the indice of the column that contains the contrast_logFC data

    if (length(indc) == 0) {
      message(
        "Please rename the logFC column to include the contrast evaluated."
      )
    } else {
      old_contrast <- colnames(diff_dat)[indc]
    }
    # actually flip contrast
    if (flip_vplot == TRUE) {
      # get the indice of the contrast to flip
      indcc <- match(old_contrast, colnames(diff_dat))
      # create flipped contrast label
      splt1 <- strsplit(old_contrast, "_") # split by underline symbol to isolate the contrast name
      splt2 <- strsplit(splt1[[1]][1], "-") # split the contrast name in the respective components
      flipped_contrast <- paste(splt2[[1]][2], splt2[[1]][1], sep = "-") # flip contrast name
      new_contrast_label <- paste(flipped_contrast, c("logFC"), sep = "_")
      # rename contrast column to the flipped contrast
      colnames(diff_dat)[indcc] <- new_contrast_label
      # flip the contrast data around y-axis
      diff_dat[, indcc] <- -diff_dat[indcc]
    } else {
      new_contrast_label <- old_contrast
    }
    plot_change_colnames <- c(plot_change_colnames, new_contrast_label)
    plot_signif_colnames <- c(plot_signif_colnames, pvalcol)

    filtered_features <- feature_ids[
      diff_dat[, pvalcol] < signif_threshold &
        abs(diff_dat[, new_contrast_label]) > change_threshold
    ]
    repeated_column <- rep(contrast, length(filtered_features))

    ## If param empty or FALSE, fill it with default value.
    if (
      is.null(add_deg_columns) ||
        length(add_deg_columns) == 0 ||
        isFALSE(add_deg_columns)
    ) {
      add_deg_columns <- c("FC", "logFC", "tstat", "pval", "adjpval")
    }

    if (all(add_deg_columns == "none")) {
      new_df <- data.frame(filtered_features, repeated_column)
      names(new_df) <- c(feature_id_colname, "Contrast")
    } else {
      add_deg_columns <- setdiff(add_deg_columns, "none")
      out_columns <- paste(contrast, add_deg_columns, sep = "_")
      deg <- diff_dat[, c(feature_id_colname, out_columns)]
      names(deg)[1] <- feature_id_colname
      new_df <- data.frame(filtered_features, repeated_column) |>
        dplyr::left_join(deg, by = c("filtered_features" = feature_id_colname))
      names(new_df) <- c(feature_id_colname, "Contrast", add_deg_columns)
    }

    df_out1 <- new_df
    df_outs[[contrast]] <- df_out1
    ### PH: END Build table for Volcano plot
  }

  # Only pass custom labels when the summary label toggles request them.
  plot_additional_labels <- ""
  if (isTRUE(add_features) || isTRUE(label_features)) {
    plot_additional_labels <- custom_gene_list
  }
  plot_titles <- gsub("_logFC$", "", plot_change_colnames)
  plot_titles <- if (identical(title, "Volcano Plots")) plot_titles else title

  message("\nRunning Enhanced Volcano:")
  volcano_enhanced_result <- plot_volcano_enhanced(
    diff_dat,
    feature_id_colname = feature_id_colname,
    signif_colname = plot_signif_colnames,
    signif_threshold = signif_threshold,
    change_colname = plot_change_colnames,
    change_threshold = change_threshold,
    value_to_sort_the_output_dataset = value_to_sort_the_output_dataset,
    num_features_to_label = num_features_to_label,
    label_features = label_features,
    custom_gene_list = plot_additional_labels,
    label_significant_features_only = label_significant_features_only,
    label_font_size = label_font_size,
    draw_connectors = draw_connectors,
    change_sig_name = change_sig_name,
    change_lfc_name = change_lfc_name,
    title = plot_titles,
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
    image_width = image_width * dpi,
    image_height = image_height * dpi,
    dpi = dpi,
    use_default_grid_layout = use_default_grid_layout,
    number_of_rows_in_grid_layout = number_of_rows_in_grid_layout,
    scale_image_to_grid = FALSE,
    print_plots = FALSE,
    save_plots = FALSE,
    plots_subdir = plots_subdir,
    plot_filename = plot_filename
  )

  plots_list <- attr(volcano_enhanced_result, "plots")
  if (length(plots_list) > 0) {
    nplots <- length(plots_list)
    if (nplots > 1) {
      abort_packages_not_installed("patchwork")
      if (
        isTRUE(use_default_grid_layout) ||
          is.null(number_of_rows_in_grid_layout)
      ) {
        nrows <- ceiling(nplots / ceiling(sqrt(nplots)))
      } else {
        nrows <- number_of_rows_in_grid_layout
      }
      if (
        !is.numeric(nrows) || length(nrows) != 1 || is.na(nrows) || nrows < 1
      ) {
        nrows <- 1
      }
      nrows <- as.integer(nrows)
      ncols <- ceiling(nplots / nrows)
      plot_obj <- patchwork::wrap_plots(plots_list, nrow = nrows)
    } else {
      nrows <- 1
      ncols <- 1
      plot_obj <- plots_list[[1]]
    }
    print_or_save_plot(
      plot_obj,
      filename = file.path(plots_subdir, plot_filename),
      print_plots = print_plots,
      save_plots = save_plots,
      units = "px",
      width = image_width * dpi * ncols,
      height = image_height * dpi * nrows,
      dpi = dpi,
      device = graphics_device
    )
  }

  df_out <- unique(do.call("rbind", df_outs))

  return(df_out)
}
