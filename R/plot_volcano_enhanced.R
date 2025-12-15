#' @rdname plot_volcano_enhanced
#' @name plot_volcano_enhanced
#' @export
plot_volcano_enhanced <- S7::new_generic("plot_volcano_enhanced", "moo_diff", function(moo_diff,
                                                                                       feature_id_colname = NULL,
                                                                                       signif_colname = c("B-A_adjpval", "B-C_adjpval"),
                                                                                       signif_threshold = 0.05,
                                                                                       change_colname = c("B-A_logFC", "B-C_logFC"),
                                                                                       change_threshold = 1.0,
                                                                                       value_to_sort_the_output_dataset = "p-value",
                                                                                       num_features_to_label = 30,
                                                                                       use_only_addition_labels = FALSE,
                                                                                       additional_labels = "",
                                                                                       is_red = TRUE,
                                                                                       lab_size = 4,
                                                                                       change_sig_name = "p-value",
                                                                                       change_lfc_name = "log2FC",
                                                                                       title = "Volcano Plots",
                                                                                       use_custom_lab = FALSE,
                                                                                       ylim = 0,
                                                                                       custom_xlim = "",
                                                                                       xlim_additional = 0,
                                                                                       ylim_additional = 0,
                                                                                       axis_lab_size = 24,
                                                                                       point_size = 2,
                                                                                       image_width = 3000,
                                                                                       image_height = 3000,
                                                                                       dpi = 300,
                                                                                       interactive_plots = FALSE,
                                                                                       print_plots = options::opt("print_plots"),
                                                                                       save_plots = options::opt("save_plots"),
                                                                                       plots_subdir = "diff",
                                                                                       plot_filename = "volcano_enhanced.png") {
  return(S7::S7_dispatch())
})

#' @rdname plot_volcano_enhanced
#' @name plot_volcano_enhanced
#' @export
S7::method(plot_volcano_enhanced, multiOmicDataSet) <- function(moo_diff,
                                                                feature_id_colname = NULL,
                                                                signif_colname = c("B-A_adjpval", "B-C_adjpval"),
                                                                signif_threshold = 0.05,
                                                                change_colname = c("B-A_logFC", "B-C_logFC"),
                                                                change_threshold = 1.0,
                                                                value_to_sort_the_output_dataset = "p-value",
                                                                num_features_to_label = 30,
                                                                use_only_addition_labels = FALSE,
                                                                additional_labels = "",
                                                                is_red = TRUE,
                                                                lab_size = 4,
                                                                change_sig_name = "p-value",
                                                                change_lfc_name = "log2FC",
                                                                title = "Volcano Plots",
                                                                use_custom_lab = FALSE,
                                                                ylim = 0,
                                                                custom_xlim = "",
                                                                xlim_additional = 0,
                                                                ylim_additional = 0,
                                                                axis_lab_size = 24,
                                                                point_size = 2,
                                                                image_width = 3000,
                                                                image_height = 3000,
                                                                dpi = 300,
                                                                interactive_plots = FALSE,
                                                                print_plots = options::opt("print_plots"),
                                                                save_plots = options::opt("save_plots"),
                                                                plots_subdir = "diff",
                                                                plot_filename = "volcano_enhanced.png") {
  return(
    join_dfs_wide() %>%
      plot_volcano_enhanced(
        feature_id_colname,
        signif_colname,
        signif_threshold,
        change_colname,
        change_threshold,
        value_to_sort_the_output_dataset,
        num_features_to_label,
        use_only_addition_labels,
        additional_labels,
        is_red,
        lab_size,
        change_sig_name,
        change_lfc_name,
        title,
        use_custom_lab,
        ylim,
        custom_xlim,
        xlim_additional,
        ylim_additional,
        axis_lab_size,
        point_size,
        image_width,
        image_height,
        dpi,
        interactive_plots,
        print_plots,
        save_plots,
        plots_subdir,
        plot_filename,
      )
  )
}

#' Enhanced Volcano Plot
#'
#' Uses [Bioconductor's Enhanced Volcano
#' Plot](https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html).
#'
#' @inheritParams option_params
#' @inheritParams filter_counts
#'
#' @param moo_diff Differential expression analysis result from one or more contrasts. This must be a dataframe.
#' @param signif_colname column name of significance values (e.g., adjusted p-values or FDR). This column will be used
#'   to determine which points are considered significant in the volcano plot.
#' @param signif_threshold Numeric value specifying the significance cutoff for p-values (i.e. filters on
#'   `signif_colname`)
#' @param change_colname column name of fold change values.
#' @param change_threshold Numeric value specifying the fold change cutoff for significance (i.e. filters on
#'   `change_colname`)
#' @param value_to_sort_the_output_dataset How to sort the output dataset. Options are "fold-change" or "p-value".
#' @param num_features_to_label Number of top features/genes to label in the volcano plot. Default is 30.
#' @param use_only_addition_labels If `TRUE`, only the additional labels specified in `additional_labels` will be used
#'   for labeling in the volcano plot, ignoring the top features.
#' @param additional_labels comma-separated string of feature names or IDs to include in the volcano plot.
#' @param is_red Logical. If TRUE, highlights points in red.
#' @param lab_size Size of the labels in the volcano plot.
#' @param change_sig_name Name for the significance column in the plot. Default is "p-value".
#' @param change_lfc_name Name for the fold change column in the plot. Default is "log2FC".
#' @param title Title of the plot. Default is "Volcano Plots".
#' @param use_custom_lab If TRUE, uses custom labels for the plot (set by `change_sig_name` and `change_lfc_name`)
#' @param ylim Y-axis limits for the plot.
#' @param custom_xlim Custom X-axis limits for the plot.
#' @param xlim_additional Additional space to add to the X-axis limits.
#' @param ylim_additional Additional space to add to the Y-axis limits.
#' @param axis_lab_size Size of the axis labels.
#' @param point_size Size of the points in the plot.
#' @param image_width output image width in pixels - only used if save_plots is TRUE
#' @param image_height output image height in pixels - only used if save_plots is TRUE
#' @param dpi dots-per-inch of the output image (see `ggsave()`) - only used if save_plots is TRUE
#' @param plot_filename plot output filename - only used if save_plots is TRUE
#'
#' @export
#' @keywords plotters volcano
#'
#' @examples
#' plot_volcano_enhanced(nidap_deg_analysis, print_plots = TRUE)
#'
#' @rdname plot_volcano_enhanced
#' @name plot_volcano_enhanced
S7::method(plot_volcano_enhanced, S7::class_data.frame) <- function(moo_diff,
                                                                    feature_id_colname = NULL,
                                                                    signif_colname = c("B-A_adjpval", "B-C_adjpval"),
                                                                    signif_threshold = 0.05,
                                                                    change_colname = c("B-A_logFC", "B-C_logFC"),
                                                                    change_threshold = 1.0,
                                                                    value_to_sort_the_output_dataset = "p-value",
                                                                    num_features_to_label = 30,
                                                                    use_only_addition_labels = FALSE,
                                                                    additional_labels = "",
                                                                    is_red = TRUE,
                                                                    lab_size = 4,
                                                                    change_sig_name = "p-value",
                                                                    change_lfc_name = "log2FC",
                                                                    title = "Volcano Plots",
                                                                    use_custom_lab = FALSE,
                                                                    ylim = 0,
                                                                    custom_xlim = "",
                                                                    xlim_additional = 0,
                                                                    ylim_additional = 0,
                                                                    axis_lab_size = 24,
                                                                    point_size = 2,
                                                                    image_width = 3000,
                                                                    image_height = 3000,
                                                                    dpi = 300,
                                                                    interactive_plots = FALSE,
                                                                    print_plots = options::opt("print_plots"),
                                                                    save_plots = options::opt("save_plots"),
                                                                    plots_subdir = "diff",
                                                                    plot_filename = "volcano_enhanced.png") {
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

  rank <- list()
  plots_list <- list()

  # user can select multiple comparisons to create volcano plots
  for (i in seq_along(change_colname)) {
    ### PH: START Build table for Volcano plot

    lfccol <- change_colname[i]
    sigcol <- signif_colname[i]
    columns_of_interest <- c(label_col, change_colname[i], signif_colname[i])
    df <- diff_dat %>%
      dplyr::select(tidyselect::one_of(columns_of_interest)) %>%
      dplyr::mutate(!!rlang::sym(lfccol) := tidyr::replace_na(!!rlang::sym(lfccol), 0)) %>%
      dplyr::mutate(!!rlang::sym(sigcol) := tidyr::replace_na(!!rlang::sym(sigcol), 1))
    # mutate(.data[[lfc.col[i]]] = replace_na(.data[[lfc.col[i]]], 0)) %>%
    # mutate(.data[[sig.col[i]]] = replace_na(.data[[sig.col[i]]], 1))
    if (use_custom_lab == TRUE) {
      if (nchar(change_lfc_name) == 0) {
        lfc_name <- change_colname[i]
      }
      if (nchar(change_sig_name) == 0) {
        sig_name <- signif_colname[i]
      }
      colnames(df) <- c(label_col, change_lfc_name, sig_name)
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
    if (value_to_sort_the_output_dataset == "fold-change") {
      df <- df %>% dplyr::arrange(dplyr::desc(.data[[lfc_name]]))
    } else if (value_to_sort_the_output_dataset == "p-value") {
      df <- df %>% dplyr::arrange(.data[[sig_name]])
    }

    if (is_red) {
      df_sub <- df[df[[sigcol]] <= signif_threshold & abs(df[[lfccol]]) >= change_threshold, ]
    } else {
      df_sub <- df
    }

    genes_to_label <- as.character(df_sub[1:num_features_to_label, label_col])
    #        additional_labels <- unlist(str_split(additional_labels,","))
    ## Modifying Additional Labels List:
    ## Replace commas with spaces and split the string
    split_values <- unlist(strsplit(gsub(",", " ", additional_labels), " "))
    additional_labels <- split_values[split_values != ""]

    filter <- additional_labels %in% df[, label_col]
    missing_labels <- additional_labels[!filter]
    additional_labels <- additional_labels[filter]

    if (length(missing_labels) > 0) {
      message(glue::glue((
        "Could not find missing labels:\t{paste(missing_labels, collapse = ', ')}"
      )))
    }

    if (use_only_addition_labels) {
      genes_to_label <- additional_labels
    } else {
      genes_to_label <- unique(append(genes_to_label, additional_labels))
    }

    significant <- vector(length = nrow(df))
    significant[] <- "Not significant"
    significant[which(abs(df[, 2]) > change_threshold)] <- "Fold change only"
    significant[which(df[, 3] < signif_threshold)] <- "Significant only"
    significant[which(abs(df[, 2]) > change_threshold & df[, 3] < signif_threshold)] <- "Significant and fold change"

    ### PH: END Build table for Volcano plot


    ### PH: START Create Volcano plot

    ### PH: Set Axis limits - Unique feature to this plot that should be included with any Volcano plot function
    ##############################

    ## Y-axis range change:
    # fix pvalue == 0
    shapeCustom <- rep(19, nrow(df))
    maxy <- max(-log10(df[[sig_name]]), na.rm = TRUE)
    if (ylim > 0) {
      maxy <- ylim
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
    if (custom_xlim == "") {
      xlim <- c(
        floor(min(df[, lfc_name])) - xlim_additional,
        ceiling(max(df[, lfc_name])) + xlim_additional
      )
    } else if (grepl(",", custom_xlim) == FALSE) {
      xlim <- c(-1 * as.numeric(trimws(custom_xlim)), as.numeric(trimws(custom_xlim)))
    } else {
      split_values <- strsplit(custom_xlim, ",")[[1]]

      # Trim whitespace and convert to numeric values
      x_min <- as.numeric(trimws(split_values[1]))
      x_max <- as.numeric(trimws(split_values[2]))

      xlim <- c(x_min, x_max)
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
    } else {
      ylab <- bquote(~ -Log[10] ~ "p-value")
    }
    if (use_custom_lab) {
      if (lfc_name != change_colname[i]) {
        xlab <- gsub("_", " ", lfc_name)
      }
      if (sig_name != signif_colname[i]) {
        ylab <- gsub("_", " ", sig_name)
      }
    }


    volcano_plot <- EnhancedVolcano::EnhancedVolcano(
      df,
      x = lfc_name,
      y = sig_name,
      lab = df[, label_col],
      selectLab = genes_to_label,
      title = title,
      # CHANGE NW: See line 78
      subtitle = group,
      xlab = xlab,
      ylab = ylab,
      xlim = xlim,
      ylim = c(0, maxy + ylim_additional),
      pCutoff = signif_threshold,
      FCcutoff = change_threshold,
      axisLabSize = axis_lab_size,
      labSize = lab_size,
      pointSize = point_size,
      shapeCustom = shapeCustom
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
        title = title,
        # CHANGE NW: See line 78
        subtitle = group,
        xlab = xlab,
        ylab = ylab,
        xlim = xlim,
        ylim = c(0, maxy + ylim_additional),
        pCutoff = signif_threshold,
        FCcutoff = change_threshold,
        axisLabSize = axis_lab_size,
        labSize = lab_size,
        pointSize = point_size,
        shapeCustom = shapeCustom
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
  plot_patchwork <- patchwork::wrap_plots(plots_list)
  print_or_save_plot(
    plot_patchwork,
    filename = file.path(plots_subdir, plot_filename),
    print_plots = print_plots,
    save_plots = save_plots,
    units = "px",
    width = image_width,
    height = image_height,
    dpi = dpi
  )

  df_final <- cbind(diff_dat, do.call(cbind, rank))
  return(df_final)
}
