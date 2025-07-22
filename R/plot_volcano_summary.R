#' Volcano Plot - Summary [CC
#' Produces one volcano plot for each tested contrast in the input DEG table.
#'
#' It can be sorted by either fold change, t-statistic, or p-value. The returned
#' dataset includes one row for each significant gene in each contrast, and
#' contains columns from the DEG analysis of that contrast as well as columns
#' useful to the Venn diagram template downstream.
#'
#' @inheritParams option_params
#' @inheritParams plot_volcano_enhanced
#' @inheritParams filter_counts
#'
#' @param add_features Add custom_gene_list To Labels. Set TRUE when you want to label a specific set of features (features) in the "custom_gene_list" parameter" IN ADDITION to the number of features you set in the "Number of Features to Label" parameter.
#' @param label_features Select TRUE when you want to label ONLY a specific list of features(features) given in the "custom_gene_list" parameter.
#' @param custom_gene_list Provide a list of features (comma separated) to be labeled on the volcano plot. You must toggle one of the following ON to see these labels: "Add features" or "Label Only My Feature List".
#' @param default_label_color Set the color for the text used to add feature (gene) name labels to points.
#' @param custom_label_color Set the color for the specific list of features (features) provided in the "Feature List" parameter.
#' @param label_x_adj adjust position of the labels on the x-axis. Default: 0.2
#' @param label_y_adj adjust position of the labels on the y-axis. Default: 0.2
#' @param line_thickness Set the thickness of the lines in the plot. Default: 0.5
#' @param label_font_size Set the font size of the labels. Default: 4
#' @param label_font_type Set the font type of the labels. Default: 1
#' @param displace_feature_labels Set to TRUE to displace gene labels. Default: FALSE. Set TRUE if you want to displace the feature (gene) label for a specific set of features. Make sure to use custom x- and y- limits and give sufficient space for displacement; otherwise other labels than the desired ones will appear displaced.
#' @param custom_gene_list_special_label_displacement Provide a list of features (comma separated) for which you want special displacement of the feature label.
#' @param special_label_displacement_x_axis Displacement of the feature label on the x-axis. Default: 2
#' @param special_label_displacement_y_axis Displacement of the feature label on the y-axis. Default: 2
#' @param color_of_signif_threshold_line Color of the significance threshold line. Default: "blue"
#' @param color_of_non_significant_features Color of the non-significant features. Default: "black"
#' @param color_of_logfold_change_threshold_line Color of the log fold change threshold line. Default: "red"
#' @param color_of_features_meeting_only_signif_threshold Color of the features that meet only the significance threshold. Default: "lightgoldenrod2"
#' @param color_for_features_meeting_pvalue_and_foldchange_thresholds Color of the features that meet both the p-value and fold change thresholds. Default: "red"
#' @param flip_vplot Set to TRUE to flip the fold change values so that the volcano plot looks like a comparison was B-A. Default: FALSE
#' @param use_default_x_axis_limit Set to TRUE to use the default x-axis limit. Default: TRUE
#' @param x_axis_limit Custom x-axis limit. Default: c(-5, 5)
#' @param use_default_y_axis_limit Set to TRUE to use the default y-axis limit. Default: TRUE
#' @param y_axis_limit Custom y-axis limit. Default: c(0, 10)
#' @param point_size Size of the points in the plot. Default: 1
#' @param add_deg_columns Add additional columns from the DEG analysis to the output dataset. Default: FALSE
#' @param use_default_grid_layout Set to TRUE to use the default grid layout. Default: TRUE
#' @param number_of_rows_in_grid_layout Number of rows in the grid layout. Default: 1
#' @param aspect_ratio Aspect ratio of the output image. Default: 4/3
#' @param graphics_device passed to `ggsave(device)`. Default: `grDevices::png`
#' @param plot_filename Filename for the output plot. Default: "volcano_plot.png"
#'
#' @export
#' @keywords plotters volcano
#'
#' @examples
#' plot_volcano_summary(nidap_deg_analysis, print_plots = TRUE)
#'
plot_volcano_summary <- function(moo_diff,
                                 feature_id_colname = NULL,
                                 signif_colname = "pval",
                                 signif_threshold = 0.05,
                                 change_threshold = 1,
                                 value_to_sort_the_output_dataset = "t-statistic",
                                 num_features_to_label = 30,
                                 add_features = FALSE,
                                 label_features = FALSE,
                                 custom_gene_list = "",
                                 default_label_color = "black",
                                 custom_label_color = "green3",
                                 label_x_adj = 0.2,
                                 label_y_adj = 0.2,
                                 line_thickness = 0.5,
                                 label_font_size = 4,
                                 label_font_type = 1,
                                 displace_feature_labels = FALSE,
                                 custom_gene_list_special_label_displacement = "",
                                 special_label_displacement_x_axis = 2,
                                 special_label_displacement_y_axis = 2,
                                 color_of_signif_threshold_line = "blue",
                                 color_of_non_significant_features = "black",
                                 color_of_logfold_change_threshold_line = "red",
                                 color_of_features_meeting_only_signif_threshold = "lightgoldenrod2",
                                 color_for_features_meeting_pvalue_and_foldchange_thresholds = "red",
                                 flip_vplot = FALSE,
                                 use_default_x_axis_limit = TRUE,
                                 x_axis_limit = 5,
                                 use_default_y_axis_limit = TRUE,
                                 y_axis_limit = 10,
                                 point_size = 2,
                                 add_deg_columns = c("FC", "logFC", "tstat", "pval", "adjpval"),
                                 graphics_device = grDevices::png,
                                 image_width = 15,
                                 image_height = 15,
                                 dpi = 300,
                                 use_default_grid_layout = TRUE,
                                 number_of_rows_in_grid_layout = 1,
                                 aspect_ratio = 0,
                                 plot_filename = "volcano.png",
                                 print_plots = options::opt("print_plots"),
                                 save_plots = options::opt("save_plots"),
                                 plots_subdir = "figures") {
  abort_packages_not_installed("patchwork")
  diff_dat <- as.data.frame(moo_diff)

  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##
  ### PH
  # Input - DEG table from Limma DEG template
  # Output - Venn Diagrams for selected Comparisons + Simplified DEG table for selected Comparisons (Only used for Venn Diagram)
  # Purpose - Create Multiple Venn Diagrams
  ## Can we use Visualizations from Advnaced Volcano function used in the stand alone Volcano plot?

  # Basic Parameters:
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(diff_dat)[1]
  }

  #  Identify all contrasts in DEG output table
  volcols <- colnames(diff_dat)
  statcols <- volcols[grepl("logFC", volcols)]
  contrasts <- unique(gsub("_logFC", "", statcols))

  Plots <- list()
  df_outs <- list()

  #  Create Volcano for each DEG comparison
  for (contrast in contrasts) {
    ### PH: START Build table for Volcano plot
    message(paste0("Preparing table for contrast: ", contrast))
    lfccol <- paste0(contrast, "_logFC")
    pvalcol <- paste0(contrast, "_", signif_colname)
    tstatcol <- paste0(contrast, "_", "tstat")

    message(paste0("Fold change column: ", lfccol))
    message(paste0(signif_colname, " column: ", pvalcol))

    if (value_to_sort_the_output_dataset == "fold-change") {
      diff_dat %>% dplyr::arrange(dplyr::desc(abs(diff_dat[, lfccol]))) -> diff_dat
    } else if (value_to_sort_the_output_dataset == "p-value") {
      diff_dat %>% dplyr::arrange(diff_dat[, pvalcol]) -> diff_dat
    } else if (value_to_sort_the_output_dataset == "t-statistic") {
      diff_dat %>% dplyr::arrange(dplyr::desc(abs(diff_dat[, tstatcol]))) -> diff_dat
    }

    ## optional Parameter: Provide a list of features to label on Volcano plot
    ## work with a list of features
    if (add_features == T) {
      gl <- trimws(unlist(strsplit(c(custom_gene_list), ",")), which = c("both"))
      ind <- match(gl, diff_dat$Gene) # get the indices of the listed features
      custom_gene_list_ind <- c(1:num_features_to_label, ind) # when list provided
      color_gene_label <- c(rep(c(default_label_color), num_features_to_label), rep(c(custom_label_color), length(ind)))
    } else if (label_features == T) {
      gl <- trimws(unlist(strsplit(c(custom_gene_list), ",")), which = c("both")) # unpack the gene list provided by the user and remove white spaces
      ind <- match(gl, diff_dat$Gene) # get the indices of the listed features
      custom_gene_list_ind <- ind # when list provided
      color_gene_label <- rep(c(custom_label_color), length(ind))
    } else {
      if (num_features_to_label > 0) {
        custom_gene_list_ind <- 1:num_features_to_label # if no list provided label the number of features given by the user
        color_gene_label <- rep(c(default_label_color), num_features_to_label)
      } else if (num_features_to_label == 0) {
        custom_gene_list_ind <- 0
      }
    }


    ## optional Parameter: IF DEG was set up A-B User can Flip FC values so that Volcano plot looks like comparison was B-A
    ## flip contrast section
    indc <- which(colnames(diff_dat) == lfccol) # get the indice of the column that contains the contrast_logFC data

    if (length(indc) == 0) {
      message("Please rename the logFC column to include the contrast evaluated.")
    } else {
      old_contrast <- colnames(diff_dat)[indc]
    }
    # actually flip contrast
    if (flip_vplot == T) {
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

    filtered_features <- diff_dat$Gene[diff_dat[, pvalcol] < signif_threshold &
      abs(diff_dat[, new_contrast_label]) > change_threshold]
    repeated_column <- rep(contrast, length(filtered_features))
    ## If param empty, fill it with default value.
    if (length(add_deg_columns) == 0) {
      add_deg_columns <- c("FC", "logFC", "tstat", "pval", "adjpval")
    } else if (all(add_deg_columns == "none")) {
      new_df <- data.frame(filtered_features, repeated_column)
      names(new_df) <- c(feature_id_colname, "Contrast")
    } else {
      add_deg_columns <- setdiff(add_deg_columns, "none")
      out_columns <- paste(contrast, add_deg_columns, sep = "_")
      deg <- diff_dat[, c(feature_id_colname, out_columns)]
      names(deg)[1] <- feature_id_colname
      new_df <- data.frame(filtered_features, repeated_column) %>% dplyr::left_join(deg, by = c("filtered_features" = feature_id_colname))
      names(new_df) <- c(feature_id_colname, "Contrast", add_deg_columns)
    }

    df_out1 <- new_df
    df_outs[[contrast]] <- df_out1

    ### PH: END Build table for Volcano plot


    ### PH: START Make plot - Can we use Enhanced volcano function from other template to make figure instead of ggplot shown here

    message(paste0(
      "Total number of features included in volcano plot: ",
      nrow(diff_dat)
    ))
    ## special nudge/repel of specific features
    if (displace_feature_labels) {
      gn <- trimws(unlist(strsplit(
        c(custom_gene_list_special_label_displacement), ","
      )), which = c("both"))
      ind_gn <- match(gn, diff_dat$Gene[custom_gene_list_ind]) # get the indices of the listed features
      nudge_x_all <- rep(c(0.2), length(diff_dat$Gene[custom_gene_list_ind]))
      nudge_y_all <- rep(c(0.2), length(diff_dat$Gene[custom_gene_list_ind]))
      nudge_x_all[ind_gn] <- c(special_label_displacement_x_axis)
      nudge_y_all[ind_gn] <- c(special_label_displacement_y_axis)
    } else {
      nudge_x_all <- label_x_adj
      nudge_y_all <- label_y_adj
    }



    # set plot parameters
    if (use_default_y_axis_limit) {
      negative_log10_p_values <- -log10(diff_dat[, pvalcol])
      ymax <- ceiling(max(negative_log10_p_values[is.finite(negative_log10_p_values)]))
    } else {
      ymax <- y_axis_limit
    }
    if (use_default_x_axis_limit) {
      xmax1 <- ceiling(max(diff_dat[, lfccol]))
      xmax2 <- ceiling(max(-diff_dat[, lfccol]))
      xmax <- max(xmax1, xmax2)
    } else {
      xmax <- x_axis_limit
    }


    grm <- diff_dat[, c(new_contrast_label, pvalcol)]
    grm[, "neglogpval"] <- -log10(diff_dat[, pvalcol])
    colnames(grm) <- c("FC", "pval", "neglogpval")
    # message(grm[custom_gene_list_ind, ])
    p <- ggplot2::ggplot(grm, ggplot2::aes(x = !!rlang::sym("FC"), y = !!rlang::sym("neglogpval"))) + # modified by RAS
      ggplot2::theme_classic() +
      ggplot2::geom_point(color = color_of_non_significant_features, size = point_size) +
      ggplot2::geom_vline(
        xintercept = c(-change_threshold, change_threshold),
        color = color_of_logfold_change_threshold_line,
        alpha = 1.0
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(signif_threshold),
        color = color_of_signif_threshold_line,
        alpha = 1.0
      ) +
      ggplot2::geom_point(
        data = grm[diff_dat[, pvalcol] < signif_threshold, ],
        color = color_of_features_meeting_only_signif_threshold,
        size = point_size
      ) +
      ggplot2::geom_point(
        data = grm[diff_dat[, pvalcol] < signif_threshold &
          abs(grm[, "FC"]) > change_threshold, ],
        color = color_for_features_meeting_pvalue_and_foldchange_thresholds,
        size = point_size
      ) +
      ggrepel::geom_text_repel(
        data = grm[custom_gene_list_ind, ],
        label = diff_dat$Gene[custom_gene_list_ind],
        color = color_gene_label,
        fontface = label_font_type,
        nudge_x = nudge_x_all,
        nudge_y = nudge_y_all,
        size = label_font_size,
        segment.size = line_thickness
      ) +
      ggplot2::xlim(-xmax, xmax) +
      ggplot2::ylim(0, ymax) +
      ggplot2::xlab(new_contrast_label) +
      ggplot2::ylab(pvalcol)

    if (aspect_ratio > 0) {
      p <- p + ggplot2::coord_fixed(ratio = aspect_ratio)
    }

    Plots[[contrast]] <- p
    ### PH: END Make plot - Can we use Enhanced volcano function from other template to make figure instead of ggplot shown here
  }


  ## Print plots
  nplots <- length(Plots)
  if (use_default_grid_layout) {
    nrows <- ceiling(nplots / ceiling(sqrt(nplots)))
  } else {
    nrows <- number_of_rows_in_grid_layout
  }
  plot_patchwork <- patchwork::wrap_plots(Plots, nrow = nrows)
  print_or_save_plot(
    plot_patchwork,
    filename = file.path(plots_subdir, plot_filename),
    print_plots = print_plots,
    save_plots = save_plots,
    graphics_device = graphics_device
  )

  df_out <- unique(do.call("rbind", df_outs))

  return(df_out)
}
