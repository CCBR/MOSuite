#' Normalize counts
#'
#' @inheritParams filter_counts
#' @param input_in_log_counts set this to `TRUE` if counts are already log2-transformed
#' @param voom_normalization_method Normalization method to be applied to the logCPM values when using `limma::voom`
#'
#' @return `multiOmicDataSet` with normalized counts
#' @export
#'
#' @examples
#' moo <- multiOmicDataSet(
#'   sample_meta_dat = as.data.frame(nidap_sample_metadata),
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
                             number_of_histogram_legend_columns = 6,
                             number_of_image_rows = 2,
                             colors_for_plots = c(
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
                             print_plots = FALSE,
                             interactive_plots = FALSE) {
  counts_dat <- moo@counts[[count_type]] %>% as.data.frame()
  sample_metadata <- moo@sample_meta %>% as.data.frame()

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
  if (isTRUE(print_plots)) {
    pca_plot <- plot_pca(
      counts_dat = df.voom,
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
    print(pca_plot)
    hist_plot <- plot_histogram(
      counts_dat = df.voom,
      sample_metadata = sample_metadata,
      sample_id_colname = sample_id_colname,
      feature_id_colname = feature_id_colname,
      group_colname = group_colname,
      label_colname = label_colname,
      color_values = colors_for_plots,
      x_axis_label = "Normalized Counts"
    ) + ggplot2::labs(caption = "normalized counts")
    print(hist_plot)
    corHM_plot <- plot_heatmap(
      counts_dat = df.filt,
      sample_metadata = sample_metadata,
      sample_id_colname = sample_id_colname,
      feature_id_colname = feature_id_colname,
      group_colname = group_colname,
      label_colname = label_colname,
      color_values = colors_for_plots
    ) + ggplot2::labs(caption = "normalized counts")
    print(corHM_plot)
  }

  message(paste("Sample columns:", paste(colnames(df.voom)[!colnames(df.voom) %in% feature_id_colname]), collapse = ", "))

  if (isFALSE("norm" %in% names(moo@counts))) {
    moo@counts[["norm"]] <- list()
  }
  moo@counts[["norm"]][["voom"]] <- df.voom
  return(moo)
}
