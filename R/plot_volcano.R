# Volcano Plot - Summary [CCBR] (eeced39d-ed52-4b16-9847-4282971775e6): v393
#' Produces one volcano plot for each tested contrast in the input DEG table.
#'
#' It can be sorted by either fold change, t-statistic, or p-value. The returned dataset includes one row for each significant gene in each contrast, and contains columns from the DEG analysis of that contrast as well as columns useful to the Venn diagram template downstream.
plot_volcano_summary <- function(DEGAnalysis) {
  # image: png

  ## --------- ##
  ## Libraries ##
  ## --------- ##

  library(ggplot2)
  library(dplyr)
  library(ggrepel)


  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##
  ### PH
  # Input - DEG table from Limma DEG template
  # Output - Venn Diagrams for selected Comparisons + Simplified DEG table for selected Comparisons (Only used for Venn Diagram)
  # Purpose - Create Multiple Venn Diagrams
  ## Can we use Visualizations from Advnaced Volcano function used in the stand alone Volcano plot?


  # Basic Parameters:
  genesmat <- DEGAnalysis
  GeneNames <- "Gene"
  stattype <- "pval"
  P_val <- 0.05
  Log2FC <- 1

  # Gene name label Parameters
  value_to_sort_the_output_dataset <- "t-statistic"
  no_genes_to_label <- 30
  AddMyGenes <- FALSE
  LabelMyGenes <- FALSE
  MyGeneList <- ""
  GeneLabelColor <- "black"
  MyGenesLabelColor <- "green3"
  LabelXadj <- 0.2
  LabelYadj <- 0.2
  LineThickness <- 0.5
  LabelFontSize <- 4
  LabelFontType <- 1
  DisplaceGeneLabels <- FALSE
  GeneListSpecialLabelDisplacement <- ""
  SpecialLabelDisplacementXAxis <- 2
  SpecialLabelDisplacementYAxis <- 2

  # Plot Parameters
  ColorofPValueThresholdLine <- "blue"
  ColorofNonSignificantGenes <- "black"
  ColorofLogFoldChangeThresholdLine <- "red"
  ColorofGenesMeetingOnlyPValueThreshold <- "lightgoldenrod2"
  colorforgenesmeetingpvalueandfoldchangethresholds <- "red"
  flipVplot <- FALSE
  UseDefaultXAxisLimit <- TRUE
  XAxisLimit <- 5
  UseDefaultYAxisLimit <- TRUE
  YAxisLimit <- 10
  PointSize <- 2

  # Table Parameters:
  add_deg_columns <- c("FC", "logFC", "tstat", "pval", "adjpval")

  # Image Parameters:
  ImageOutputFormat <- "png"
  Usesvglite <- FALSE
  image_width <- 15
  image_height <- 15
  image_resolution <- 300
  UseDefaultGridLayout <- TRUE
  NumberofRowsinGridLayout <- 1
  aspect_ratio <- 0


  ## --------------- ##
  ## Error Messages ##
  ## -------------- ##

  # None so far


  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##


  # Set Visualation parameters
  if (Usesvglite & ((ImageOutputFormat == "svg"))) {
    library(svglite)
    svglite::svglite(
      file = graphicsFile,
      width = image_width,
      height = image_height,
      pointsize = 1,
      bg = "white",
    )
  } else if (ImageOutputFormat == "png") {
    png(
      filename = graphicsFile,
      width = image_width,
      height = image_height,
      units = "in",
      pointsize = 4,
      bg = "white",
      res = image_resolution,
      type = "cairo"
    )
  }




  #  Identify all contrasts in DEG output table
  volcols <- colnames(genesmat)
  print(volcols)
  statcols <- volcols[grepl("logFC", volcols)]
  contrasts <- unique(gsub("_logFC", "", statcols))

  Plots <- list()
  df_outs <- list()

  #  Create Volcano for each DEG comparison
  for (contrast in contrasts) {
    ### PH: START Build table for Volcano plot
    print(paste0("Doing contrast: ", contrast))
    lfccol <- paste0(contrast, "_logFC")
    pvalcol <- paste0(contrast, "_", stattype)
    tstatcol <- paste0(contrast, "_", "tstat")

    print(paste0("Fold change column: ", lfccol))
    print(paste0(stattype, " column: ", pvalcol))


    if (value_to_sort_the_output_dataset == "fold-change") {
      genesmat %>% dplyr::arrange(desc(abs(genesmat[, lfccol]))) -> genesmat
    } else if (value_to_sort_the_output_dataset == "p-value") {
      genesmat %>% dplyr::arrange(genesmat[, pvalcol]) -> genesmat
    } else if (value_to_sort_the_output_dataset == "t-statistic") {
      genesmat %>% dplyr::arrange(desc(abs(genesmat[, tstatcol]))) -> genesmat
    }

    ## optional Parameter: Provide a list of genes to label on Volcano plot
    ## work with a list of genes
    if (AddMyGenes == T) {
      gl <- trimws(unlist(strsplit(c(MyGeneList), ",")), which = c("both"))
      ind <- match(gl, genesmat$Gene) # get the indices of the listed genes
      gene_list_ind <- c(1:no_genes_to_label, ind) # when list provided
      color_gene_label <- c(rep(c(GeneLabelColor), no_genes_to_label), rep(c(MyGenesLabelColor), length(ind)))
    } else if (LabelMyGenes == T) {
      gl <- trimws(unlist(strsplit(c(MyGeneList), ",")), which = c("both")) # unpack the gene list provided by the user and remove white spaces
      ind <- match(gl, genesmat$Gene) # get the indices of the listed genes
      gene_list_ind <- ind # when list provided
      color_gene_label <- rep(c(MyGenesLabelColor), length(ind))
    } else {
      if (no_genes_to_label > 0) {
        gene_list_ind <- 1:no_genes_to_label # if no list provided label the number of genes given by the user
        color_gene_label <- rep(c(GeneLabelColor), no_genes_to_label)
      } else if (no_genes_to_label == 0) {
        gene_list_ind <- 0
      }
    }


    ## optional Parameter: IF DEG was set up A-B User can Flip FC values so that Volcano plot looks like comparison was B-A
    ## flip contrast section
    indc <- which(colnames(genesmat) == lfccol) # get the indice of the column that contains the contrast_logFC data


    if (length(indc) == 0) {
      print("Please rename the logFC column to include the contrast evaluated.")
    } else {
      old_contrast <- colnames(genesmat)[indc]
    }
    # actually flip contrast
    if (flipVplot == T) {
      # get the indice of the contrast to flip
      indcc <- match(old_contrast, colnames(genesmat))
      # create flipped contrast label
      splt1 <- strsplit(old_contrast, "_") # split by underline symbol to isolate the contrast name
      splt2 <- strsplit(splt1[[1]][1], "-") # split the contrast name in the respective components
      flipped_contrast <- paste(splt2[[1]][2], splt2[[1]][1], sep = "-") # flip contrast name
      new_contrast_label <- paste(flipped_contrast, c("logFC"), sep = "_")
      # rename contrast column to the flipped contrast
      colnames(genesmat)[indcc] <- new_contrast_label
      # flip the contrast data around y-axis
      genesmat[, indcc] <- -genesmat[indcc]
    } else {
      new_contrast_label <- old_contrast
    }




    filtered_genes <- genesmat$Gene[genesmat[, pvalcol] < P_val & abs(genesmat[, new_contrast_label]) > Log2FC]
    # print(filtered_genes)
    repeated_column <- rep(contrast, length(filtered_genes))
    ## If param empty, fill it with default value.
    if (length(add_deg_columns) == 0) {
      add_deg_columns <- c("FC", "logFC", "tstat", "pval", "adjpval")
    } else if (add_deg_columns == "none") {
      new_df <- data.frame(filtered_genes, repeated_column)
      names(new_df) <- c("Gene", "Contrast")
    } else {
      add_deg_columns <- setdiff(add_deg_columns, "none")
      out_columns <- paste(contrast, add_deg_columns, sep = "_")
      deg <- genesmat[, c("Gene", out_columns)]
      names(deg)[1] <- "Gene"
      new_df <- data.frame(filtered_genes, repeated_column) %>% dplyr::left_join(deg, by = c("filtered_genes" = "Gene"))
      names(new_df) <- c("Gene", "Contrast", add_deg_columns)
    }

    df_out1 <- new_df
    df_outs[[contrast]] <- df_out1

    ### PH: END Build table for Volcano plot


    ### PH: START Make plot - Can we use Enhanced volcano function from other template to make figure instead of ggplot shown here

    print(paste0("Total number of genes included in volcano plot: ", nrow(genesmat)))
    ## special nudge/repel of specific genes
    if (DisplaceGeneLabels) {
      gn <- trimws(unlist(strsplit(c(GeneListSpecialLabelDisplacement), ",")), which = c("both"))
      ind_gn <- match(gn, genesmat$Gene[gene_list_ind]) # get the indices of the listed genes
      nudge_x_all <- rep(c(0.2), length(genesmat$Gene[gene_list_ind]))
      nudge_y_all <- rep(c(0.2), length(genesmat$Gene[gene_list_ind]))
      nudge_x_all[ind_gn] <- c(SpecialLabelDisplacementXAxis)
      nudge_y_all[ind_gn] <- c(SpecialLabelDisplacementYAxis)
    } else {
      nudge_x_all <- LabelXadj
      nudge_y_all <- LabelYadj
    }



    # set plot parameters
    if (UseDefaultYAxisLimit) {
      negative_log10_p_values <- -log10(genesmat[, pvalcol])
      ymax <- ceiling(max(negative_log10_p_values[is.finite(negative_log10_p_values)]))
    } else {
      ymax <- YAxisLimit
    }
    if (UseDefaultXAxisLimit) {
      xmax1 <- ceiling(max(genesmat[, lfccol]))
      xmax2 <- ceiling(max(-genesmat[, lfccol]))
      xmax <- max(xmax1, xmax2)
    } else {
      xmax <- XAxisLimit
    }


    grm <- genesmat[, c(new_contrast_label, pvalcol)]
    grm[, "neglogpval"] <- -log10(genesmat[, pvalcol])
    colnames(grm) <- c("FC", "pval", "neglogpval")
    print(grm[gene_list_ind, ])
    p <- ggplot(
      grm,
      aes_string(x = "FC", y = "neglogpval")
    ) + # modified by RAS
      theme_classic() +
      geom_point(
        color = ColorofNonSignificantGenes,
        size = PointSize
      ) +
      geom_vline(xintercept = c(-Log2FC, Log2FC), color = ColorofLogFoldChangeThresholdLine, alpha = 1.0) +
      geom_hline(yintercept = -log10(P_val), color = ColorofPValueThresholdLine, alpha = 1.0) +
      geom_point(
        data = grm[genesmat[, pvalcol] < P_val, ],
        color = ColorofGenesMeetingOnlyPValueThreshold,
        size = PointSize
      ) +
      geom_point(
        data = grm[genesmat[, pvalcol] < P_val & abs(grm[, "FC"]) > Log2FC, ],
        color = colorforgenesmeetingpvalueandfoldchangethresholds,
        size = PointSize
      ) +
      geom_text_repel(
        data = grm[gene_list_ind, ],
        label = genesmat$Gene[gene_list_ind],
        color = color_gene_label,
        fontface = LabelFontType,
        nudge_x = nudge_x_all,
        nudge_y = nudge_y_all,
        size = LabelFontSize,
        segment.size = LineThickness
      ) +
      xlim(-xmax, xmax) +
      ylim(0, ymax) +
      xlab(new_contrast_label) +
      ylab(pvalcol)

    if (aspect_ratio > 0) {
      p <- p + coord_fixed(ratio = aspect_ratio)
    }

    # print(p)
    Plots[[contrast]] <- p
    ### PH: END Make plot - Can we use Enhanced volcano function from other template to make figure instead of ggplot shown here
  }


  ## Print plots
  Use_default_grid_layout <- UseDefaultGridLayout
  require(gridExtra)
  nplots <- length(Plots)
  if (Use_default_grid_layout) {
    nrows <- ceiling(nplots / ceiling(sqrt(nplots)))
  } else {
    nrows <- NumberofRowsinGridLayout
  }

  do.call("grid.arrange", c(Plots, nrow = nrows))
  print("done plotting")

  df_out <- unique(do.call("rbind", df_outs))
  print(head(df_out))
  print(colnames(df_out))


  return(df_out)
}


# Volcano Plot - Enhanced [CCBR] [scRNA-seq] [Bulk] (0c91aa57-0f76-4513-a063-5f9263d65727): v85
#' Implementation of Bioconductor's Enhanced Volcano Plot (v1.6.0, https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html).
plot_volcano_enhanced <- function(DEGAnalysis) {
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
  df.orig <- DEGAnalysis
  label.col <- "Gene"
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
  ## Error Messages ##
  ## -------------- ##

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
            "<br>Log2FC:", df[[lfc_name]],
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
