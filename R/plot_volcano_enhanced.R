# Volcano Plot - Enhanced [CCBR] [scRNA-seq] [Bulk] (0c91aa57-0f76-4513-a063-5f9263d65727): v85
#' Implementation of Bioconductor's Enhanced Volcano Plot (v1.6.0, https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html).
plot_volcano_enhanced <- function(DEGAnalysis, feature_id_colname = NULL) {
  # image: png

  # Changelog
  # 2022-09-14 Rearranged structure and description
  # 2020-10-29 Add support for pval == 0

  ## --------- ##
  ## Libraries ##
  ## --------- ##

  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)

  library(EnhancedVolcano)
  # For interactive plot:
  library(plotly)
  library(grid)


  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##
  ### PH
  # Input - DEG table from Limma DEG template
  # Output - Volcano plot + interactive Volcano Plot
  # Purpouse Create detailed Volcano for single Comparison

  # Basic Parameters:
  df.orig <- as.data.frame(DEGAnalysis)
  label.col <- feature_id_colname
  sig.col <- c("B-A_adjpval", "B-C_adjpval")
  pCutoff <- 0.05
  lfc.col <- c("B-A_logFC", "B-C_logFC")
  FCcutoff <- 1.0


  # Label Parameters
  value_to_sort_the_output_dataset <- "p-value"
  no_genes_to_label <- 30
  use_only_addition_labels <- FALSE
  additional_labels <- ""
  is_red <- TRUE
  labSize <- 4


  # Title and Axis labels Parameters
  change_sig_name <- "p-value"
  change_lfc_name <- "log2FC"
  title <- "Volcano Plots"
  # subtitle <- ""
  use_custom_lab <- FALSE

  # Plot Parameters
  ylim <- 0
  custom_xlim <- ""
  xlim_additional <- 0
  ylim_additional <- 0
  axisLabSize <- 24
  pointSize <- 2


  # Image Parameters
  imageWidth <- 3000
  imageHeight <- 3000
  dpi <- 300

  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##

  rank <- list()

  # user can select multiple comparisons to create volcano plots
  for (i in 1:length(lfc.col)) {
    ### PH: START Build table for Volcano plot

    lfccol <- lfc.col[i]
    sigcol <- sig.col[i]
    columns_of_interest <- c(label.col, lfc.col[i], sig.col[i])
    df <- df.orig %>%
      dplyr::select(one_of(columns_of_interest)) %>%
      mutate(!!sym(lfccol) := replace_na(!!sym(lfccol), 0)) %>%
      mutate(!!sym(sigcol) := replace_na(!!sym(sigcol), 1))
    # mutate(.data[[lfc.col[i]]] = replace_na(.data[[lfc.col[i]]], 0)) %>%
    # mutate(.data[[sig.col[i]]] = replace_na(.data[[sig.col[i]]], 1))
    if (use_custom_lab == TRUE) {
      if (nchar(change_lfc_name) == 0) {
        lfc_name <- lfc.col[i]
      }
      if (nchar(change_sig_name) == 0) {
        sig_name <- sig.col[i]
      }
      colnames(df) <- c(label.col, change_lfc_name, sig_name)
    } else {
      lfc_name <- lfc.col[i]
      sig_name <- sig.col[i]
    }

    ### PH: START Creating rank based on pvalue and fold change
    ## This is unique to this template and could be useful as a generic tool to create rankes for GSEA. Recommend extracting this function
    group <- gsub("_pval|p_val_", "", sig_name)
    rank[[i]] <- -log10(df[[sig_name]]) * sign(df[[lfc_name]])
    names(rank)[i] <- paste0("C_", group, "_rank")
    ### PH: End Creating rank based on pvalue and fold change

    cat(paste0("Genes in initial dataset: ", nrow(df), "\n"))

    # Select top genes by logFC or Significance
    if (value_to_sort_the_output_dataset == "fold-change") {
      df <- df %>% dplyr::arrange(desc(.data[[lfc_name]]))
    } else if (value_to_sort_the_output_dataset == "p-value") {
      df <- df %>% dplyr::arrange(.data[[sig_name]])
    }

    if (is_red) {
      df_sub <- df[df[[sigcol]] <= pCutoff & abs(df[[lfccol]]) >= FCcutoff, ]
    } else {
      df_sub <- df
    }

    genes_to_label <- as.character(df_sub[1:no_genes_to_label, label.col])
    #        additional_labels <- unlist(str_split(additional_labels,","))
    ## Modifying Additional Labels List:
    ## Replace commas with spaces and split the string
    split_values <- unlist(strsplit(gsub(",", " ", additional_labels), " "))
    additional_labels <- split_values[split_values != ""]

    filter <- additional_labels %in% df[, label.col]
    missing_labels <- additional_labels[!filter]
    additional_labels <- additional_labels[filter]

    if (length(missing_labels) > 0) {
      cat("Could not find:\n")
      print(missing_labels)
    }

    if (use_only_addition_labels) {
      genes_to_label <- additional_labels
    } else {
      genes_to_label <- unique(append(genes_to_label, additional_labels))
    }

    significant <- vector(length = nrow(df))
    significant[] <- "Not significant"
    significant[which(abs(df[, 2]) > FCcutoff)] <- "Fold change only"
    significant[which(df[, 3] < pCutoff)] <- "Significant only"
    significant[which(abs(df[, 2]) > FCcutoff & df[, 3] < pCutoff)] <- "Significant and fold change"
    print(table(significant))

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

    cat(paste0("Maxy: ", maxy, "\n"))
    if (maxy == Inf) {
      # Sometimes, pvalues == 0
      keep <- df[[sig_name]] > 0
      df[[sig_name]][!keep] <- min(df[[sig_name]][keep])
      shapeCustom[!keep] <- 17

      maxy <- -log10(min(df[[sig_name]][keep]))
      cat("Some p-values equal zero. Adjusting y-limits.\n")
      cat(paste0("Maxy adjusted: ", maxy, "\n"))
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
      xlim <- c(floor(min(df[, lfc_name])) - xlim_additional, ceiling(max(df[, lfc_name])) + xlim_additional)
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

    if (grepl("log", lfc.col[i])) {
      xlab <- bquote(~ Log[2] ~ "fold change")
    } else {
      xlab <- "Fold change"
    }
    if (grepl("adj", sig.col[i])) {
      ylab <- bquote(~ -Log[10] ~ "FDR")
    } else {
      ylab <- bquote(~ -Log[10] ~ "p-value")
    }
    if (use_custom_lab) {
      if (lfc_name != lfc.col[i]) {
        xlab <- gsub("_", " ", lfc_name)
      }
      if (sig_name != sig.col[i]) {
        ylab <- gsub("_", " ", sig_name)
      }
    }


    p <- EnhancedVolcano(df,
      x = lfc_name, y = sig_name,
      lab = df[, label.col],
      selectLab = genes_to_label,
      title = title, # CHANGE NW: See line 78
      subtitle <- group,
      xlab = xlab,
      ylab = ylab,
      xlim = xlim,
      ylim = c(0, maxy + ylim_additional),
      pCutoff = pCutoff,
      FCcutoff = FCcutoff,
      axisLabSize = axisLabSize,
      labSize = labSize,
      pointSize = pointSize,
      shapeCustom = shapeCustom
    )
    print(p)


    ## Creating plot that can be converted to plotly interactive plot (no labels):
    ## PH: make this feature an option not default
    p_empty <- EnhancedVolcano(df,
      x = lfc_name, y = sig_name,
      lab = rep("", nrow(df)), # Setting labels to empty strings
      selectLab = NULL,
      title = title, # CHANGE NW: See line 78
      subtitle <- group,
      xlab = xlab,
      ylab = ylab,
      xlim = xlim,
      ylim = c(0, maxy + ylim_additional),
      pCutoff = pCutoff,
      FCcutoff = FCcutoff,
      axisLabSize = axisLabSize,
      labSize = labSize,
      pointSize = pointSize,
      shapeCustom = shapeCustom
    )
    ##        print(p_empty)

    # Extract the data used for plotting
    plot_data <- ggplot_build(p_empty)$data[[1]]

    pxx <- p_empty +
      xlab("Fold Change") + # Simplify x-axis label
      ylab("Significance") + # Simplify y-axis label
      theme_minimal() +
      geom_point(
        aes(
          text = paste(
            "Gene:", df[[label.col]],
            "<br>log2fc_threshold:", df[[lfc_name]],
            "<br>P-value:", df[[sig_name]]
          ),
          colour = as.character(plot_data$colour),
          fill = as.character(plot_data$colour) # Set fill to the same as colour
        ),
        shape = 21, # Shape that supports both colour and fill
        size = 2, # Size of the points
        stroke = 0.1 # Stroke width
      ) + scale_fill_identity()

    # Add interactive hover labels for the gene names
    interactive_plot <- ggplotly(pxx, tooltip = c("text"))

    ### PH: END Create Volcano plot


    ## showing interactive plot
    grid.newpage()
    print(interactive_plot)
    grid.newpage()
    ## The end of addition
  }

  df.final <- cbind(df.orig, do.call(cbind, rank))
  return(df.final)
}

#################################################
## Global imports and functions included below ##
#################################################
