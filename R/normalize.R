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
#'   counts_lst = list(
#'     "raw" = as.data.frame(nidap_raw_counts),
#'     "clean" = as.data.frame(nidap_clean_raw_counts),
#'     "filt" = as.data.frame(nidap_filtered_counts)
#'   )
#' ) %>%
#'   normalize_counts(
#'     gene_names_column = "Gene",
#'     columns_to_include = c("Gene", "A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"),
#'     sample_names_column = "Sample",
#'     group_column = "Group",
#'     label_column = "Label"
#'   )
#' head(moo@counts[["norm"]][["voom"]])
normalize_counts <- function(moo,
                             count_type = "filt",
                             gene_names_column = "Gene",
                             columns_to_include = c("Gene", "A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"),
                             sample_names_column = "Sample",
                             group_column = "Group",
                             label_column = "Label",
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
                             colors_for_plots = c(),
                             make_plots_interactive = FALSE,
                             plot_correlation_matrix_heatmap = TRUE) {
  counts_matrix <- moo@counts[[count_type]] %>% as.data.frame()
  sample_metadata <- moo@sample_meta %>% as.data.frame()

  ### PH: START Color Picking Function
  ##################################
  ##### Color Picking Function
  ##################################

  ### Do not Functionalize This code yet.
  ### I will make a General Color picker function to create a custom Color Palette to use for all figures
  ### I think I would also like to use the metadata table to manually set colors
  colorlist <- c(
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
  )
  names(colorlist) <- c(
    "indigo",
    "carrot",
    "lipstick",
    "turquoise",
    "lavender",
    "jade",
    "coral",
    "azure",
    "green",
    "rum",
    "orange",
    "olive"
  )
  if (length(colors_for_plots) == 0) {
    colors_for_plots <- c(
      "indigo",
      "carrot",
      "lipstick",
      "turquoise",
      "lavender",
      "jade",
      "coral",
      "azure",
      "green",
      "rum",
      "orange",
      "olive"
    )
  }
  colorval <- colorlist[colors_for_plots]
  colorval <- unname(colorval) # remove names which affect ggplot
  ### PH: END Color Picking Function

  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##


  ##############################
  #### Create Generic Row names
  ##############################
  ## create unique rownames to correctly add back Annocolumns at end of template

  samples_to_include <- columns_to_include[columns_to_include %in% sample_metadata[, sample_names_column, drop = T]]
  anno_col <- columns_to_include[columns_to_include %in% sample_metadata[, sample_names_column, drop = T] == F]


  samples_to_include <- samples_to_include[samples_to_include != gene_names_column]
  samples_to_include <- samples_to_include[samples_to_include != "Gene"]
  samples_to_include <- samples_to_include[samples_to_include != "GeneName"]

  ## create unique rownames to correctly add back Annocolumns at end of template
  counts_matrix[, gene_names_column] <- paste0(counts_matrix[, gene_names_column], "_", 1:nrow(counts_matrix))


  anno_col <- c(anno_col, gene_names_column) %>% unique()
  anno_tbl <- counts_matrix[, anno_col, drop = F] %>% as.data.frame()


  df.filt <- counts_matrix[, samples_to_include]
  gene_names <- NULL
  gene_names$GeneID <- counts_matrix[, 1]

  ##############################
  #### Input Data Validation
  ##############################
  # TODO move this logic to S7 validator
  sample_metadata <- sample_metadata[match(colnames(df.filt), sample_metadata[[sample_names_column]]), ] # First match sample metadata to counts matrix
  sample_metadata <- sample_metadata[rowSums(is.na(sample_metadata)) != ncol(sample_metadata), ] # Remove empty rows
  sample_metadata <- sample_metadata[, colSums(is.na(sample_metadata)) == 0] # Remove empty columns
  rownames(sample_metadata) <- sample_metadata[[sample_names_column]]

  df.filt <- df.filt[, match(sample_metadata[[sample_names_column]], colnames(df.filt))] # Match counts matrix columns to sample metadata


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
  rownames(v$E) <- v$genes$GeneID
  as.data.frame(v$E) %>% tibble::rownames_to_column(gene_names_column) -> df.voom
  message(paste0("Total number of features included: ", nrow(df.voom)))
  ### PH: END Limma Normalization

  pca_plot <- plot_pca(v$E,
    sample_metadata,
    samples_to_include,
    samples_to_rename,
    group_column,
    label_column,
    color_values = colorval,
    principal_component_on_x_axis = principal_component_on_x_axis,
    principal_component_on_y_axis = principal_component_on_y_axis,
    legend_position_for_pca = legend_position_for_pca,
    point_size_for_pca = point_size_for_pca,
    add_label_to_pca = add_label_to_pca,
    label_font_size = label_font_size,
    label_offset_y_ = label_offset_y_,
    label_offset_x_ = label_offset_x_
  )

  histPlot <- plot_histogram(
    v$E,
    sample_metadata,
    gene_names_column,
    group_column,
    label_column,
    color_values = colorval,
    x_axis_label = "Normalized Counts"
  )
  if (plot_correlation_matrix_heatmap == TRUE) {
    corHM <- plot_heatmap(
      counts_matrix = df.filt,
      sample_metadata = sample_metadata,
      anno_colors = colorval,
      anno_column = group_column,
      label_column = label_column,
      sample_names_column = sample_names_column
    )
  }

  message("Sample columns")
  message(colnames(df.voom)[!colnames(df.voom) %in% gene_names_column])
  message("Feature Columns")
  message(colnames(anno_tbl))

  df.voom <- merge(anno_tbl, df.voom, by = gene_names_column, all.y = T)
  df.voom[, gene_names_column] <- gsub("_[0-9]+$", "", df.voom[, gene_names_column])

  if (isFALSE("norm" %in% names(moo@counts))) {
    moo@counts[["norm"]] <- list()
  }
  moo@counts[["norm"]][["voom"]] <- df.voom
  return(moo)
}
