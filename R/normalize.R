#' Normalize counts
#'
#' @inheritParams filter_counts
#' @inheritParams option_params
#'
#' @param norm_type normalization type. Default: "voom" which uses `limma::voom`.
#' @param input_in_log_counts set this to `TRUE` if counts are already log2-transformed
#' @param voom_normalization_method Normalization method to be applied to the logCPM values when using `limma::voom`
#'
#' @return `multiOmicDataSet` with normalized counts
#' @export
#'
#' @examples
#' moo <- multiOmicDataSet(
#'   sample_metadata = as.data.frame(nidap_sample_metadata),
#'   anno_dat = data.frame(),
#'   counts_lst = list(
#'     "raw" = as.data.frame(nidap_raw_counts),
#'     "clean" = as.data.frame(nidap_clean_raw_counts),
#'     "filt" = as.data.frame(nidap_filtered_counts)
#'   )
#' ) %>%
#'   normalize_counts(
#'     group_colname = "Group",
#'     label_colname = "Label"
#'   )
#' head(moo@counts[["norm"]][["voom"]])
normalize_counts <- function(moo,
                             count_type = "filt",
                             norm_type = "voom",
                             feature_id_colname = NULL,
                             samples_to_include = NULL,
                             sample_id_colname = NULL,
                             group_colname = "Group",
                             label_colname = NULL,
                             input_in_log_counts = FALSE,
                             voom_normalization_method = "quantile",
                             samples_to_rename = c(""),
                             add_label_to_pca = TRUE,
                             principal_component_on_x_axis = 1,
                             principal_component_on_y_axis = 2,
                             legend_position_for_pca = "top",
                             label_offset_x_ = 2,
                             label_offset_y_ = 2,
                             label_font_size = 3,
                             point_size_for_pca = 8,
                             color_histogram_by_group = TRUE,
                             set_min_max_for_x_axis_for_histogram = FALSE,
                             minimum_for_x_axis_for_histogram = -1,
                             maximum_for_x_axis_for_histogram = 1,
                             legend_font_size_for_histogram = 10,
                             legend_position_for_histogram = "top",
                             colors_for_plots = NULL,
                             print_plots = options::opt("print_plots"),
                             save_plots = options::opt("save_plots"),
                             interactive_plots = FALSE,
                             plots_subdir = "norm") {
  message(glue::glue("* normalizing {count_type} counts"))
  counts_dat <- moo@counts[[count_type]] %>% as.data.frame()
  sample_metadata <- moo@sample_meta %>% as.data.frame()
  plots_subdir <- file.path(plots_subdir, norm_type)
  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(counts_dat)[1]
  }
  if (is.null(samples_to_include)) {
    samples_to_include <- sample_metadata %>% dplyr::pull(sample_id_colname)
  }
  if (is.null(label_colname)) {
    label_colname <- sample_id_colname
  }
  df.filt <- counts_dat %>%
    dplyr::select(tidyselect::all_of(samples_to_include))


  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  gene_names <- NULL
  gene_names$feature_id <- counts_dat %>% dplyr::pull(feature_id_colname)

  ### PH: START Limma Normalization
  ##############################
  #### Limma Normalization
  ##############################

  # If input is in log space, linearize
  if (input_in_log_counts == TRUE) {
    x <- edgeR::DGEList(counts = 2^df.filt, genes = gene_names)
  } else {
    x <- edgeR::DGEList(counts = df.filt, genes = gene_names)
  }
  v <- limma::voom(x, normalize = voom_normalization_method)
  rownames(v$E) <- v$genes$feature_id
  df.voom <- as.data.frame(v$E) %>% tibble::rownames_to_column(feature_id_colname)
  message(paste0("Total number of features included: ", nrow(df.voom)))
  ### PH: END Limma Normalization
  if (isTRUE(print_plots) | isTRUE(save_plots)) {
    if (is.null(colors_for_plots)) {
      colors_for_plots <- moo@analyses[["colors"]][[group_colname]]
    }
    if (isTRUE(color_histogram_by_group)) {
      colors_for_histogram <- colors_for_plots
    } else {
      colors_for_histogram <- moo@analyses[["colors"]][[label_colname]]
    }
    pca_plot <- plot_pca(
      df.voom,
      sample_metadata = sample_metadata,
      sample_id_colname = sample_id_colname,
      samples_to_rename = samples_to_rename,
      group_colname = group_colname,
      label_colname = label_colname,
      color_values = colors_for_plots,
      principal_components = c(
        principal_component_on_x_axis,
        principal_component_on_y_axis
      ),
      legend_position = legend_position_for_pca,
      point_size = point_size_for_pca,
      add_label = add_label_to_pca,
      label_font_size = label_font_size,
      label_offset_y_ = label_offset_y_,
      label_offset_x_ = label_offset_x_
    ) + ggplot2::labs(caption = "normalized counts")
    hist_plot <- plot_histogram(
      df.voom,
      sample_metadata = sample_metadata,
      sample_id_colname = sample_id_colname,
      feature_id_colname = feature_id_colname,
      group_colname = group_colname,
      label_colname = label_colname,
      color_values = colors_for_histogram,
      color_by_group = color_histogram_by_group,
      x_axis_label = "Normalized Counts",
      legend_position = legend_position_for_histogram,
      legend_font_size = legend_font_size_for_histogram
    ) + ggplot2::labs(caption = "normalized counts")
    corHM_plot <- plot_corr_heatmap(
      df.filt,
      sample_metadata = sample_metadata,
      sample_id_colname = sample_id_colname,
      feature_id_colname = feature_id_colname,
      group_colname = group_colname,
      label_colname = label_colname,
      color_values = colors_for_plots
    ) + ggplot2::labs(caption = "normalized counts")

    print_or_save_plot(pca_plot,
      filename = file.path(plots_subdir, "pca.png"),
      print_plots = print_plots, save_plots = save_plots
    )
    print_or_save_plot(hist_plot,
      filename = file.path(plots_subdir, "histogram.png"),
      print_plots = print_plots, save_plots = save_plots
    )
    print_or_save_plot(corHM_plot,
      filename = file.path(plots_subdir, "corr_heatmap.png"),
      print_plots = print_plots, save_plots = save_plots
    )
  }

  message(paste("Sample columns:", paste(colnames(df.voom)[!colnames(df.voom) %in% feature_id_colname]), collapse = ", "))

  if (isFALSE("norm" %in% names(moo@counts))) {
    moo@counts[["norm"]] <- list()
  }
  moo@counts[["norm"]][[norm_type]] <- df.voom
  return(moo)
}
