#' Volcano Plot - Summary [CCBR] (eeced39d-ed52-4b16-9847-4282971775e6): v393
#' Produces one volcano plot for each tested contrast in the input DEG table.
#'
#' It can be sorted by either fold change, t-statistic, or p-value. The returned
#' dataset includes one row for each significant gene in each contrast, and
#' contains columns from the DEG analysis of that contrast as well as columns
#' useful to the Venn diagram template downstream.
S7::method(plot_volcano_summary, S7::class_data.frame) <- function(moo_counts,
                                                                   feature_id_colname = NULL, ...) {
  # image: png

  ## --------- ##
  ## Libraries ##
  ## --------- ##

  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  counts_dat <- as.data.frame(moo_counts)

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
    feature_id_colname <- colnames(counts_dat)[1]
  }

  stat_type <- "pval"
  pvalue_threshold <- 0.05
  log2fc_threshold <- 1

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
  graphicsFile <- "volcano.png"

  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##

  # Set Visualization parameters
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
  volcols <- colnames(counts_dat)
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
    pvalcol <- paste0(contrast, "_", stat_type)
    tstatcol <- paste0(contrast, "_", "tstat")

    print(paste0("Fold change column: ", lfccol))
    print(paste0(stat_type, " column: ", pvalcol))

    if (value_to_sort_the_output_dataset == "fold-change") {
      counts_dat %>% dplyr::arrange(desc(abs(counts_dat[, lfccol]))) -> counts_dat
    } else if (value_to_sort_the_output_dataset == "p-value") {
      counts_dat %>% dplyr::arrange(counts_dat[, pvalcol]) -> counts_dat
    } else if (value_to_sort_the_output_dataset == "t-statistic") {
      counts_dat %>% dplyr::arrange(desc(abs(counts_dat[, tstatcol]))) -> counts_dat
    }

    ## optional Parameter: Provide a list of genes to label on Volcano plot
    ## work with a list of genes
    if (AddMyGenes == T) {
      gl <- trimws(unlist(strsplit(c(MyGeneList), ",")), which = c("both"))
      ind <- match(gl, counts_dat$Gene) # get the indices of the listed genes
      gene_list_ind <- c(1:no_genes_to_label, ind) # when list provided
      color_gene_label <- c(rep(c(GeneLabelColor), no_genes_to_label), rep(c(MyGenesLabelColor), length(ind)))
    } else if (LabelMyGenes == T) {
      gl <- trimws(unlist(strsplit(c(MyGeneList), ",")), which = c("both")) # unpack the gene list provided by the user and remove white spaces
      ind <- match(gl, counts_dat$Gene) # get the indices of the listed genes
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
    indc <- which(colnames(counts_dat) == lfccol) # get the indice of the column that contains the contrast_logFC data

    if (length(indc) == 0) {
      print("Please rename the logFC column to include the contrast evaluated.")
    } else {
      old_contrast <- colnames(counts_dat)[indc]
    }
    # actually flip contrast
    if (flipVplot == T) {
      # get the indice of the contrast to flip
      indcc <- match(old_contrast, colnames(counts_dat))
      # create flipped contrast label
      splt1 <- strsplit(old_contrast, "_") # split by underline symbol to isolate the contrast name
      splt2 <- strsplit(splt1[[1]][1], "-") # split the contrast name in the respective components
      flipped_contrast <- paste(splt2[[1]][2], splt2[[1]][1], sep = "-") # flip contrast name
      new_contrast_label <- paste(flipped_contrast, c("logFC"), sep = "_")
      # rename contrast column to the flipped contrast
      colnames(counts_dat)[indcc] <- new_contrast_label
      # flip the contrast data around y-axis
      counts_dat[, indcc] <- -counts_dat[indcc]
    } else {
      new_contrast_label <- old_contrast
    }

    filtered_genes <- counts_dat$Gene[counts_dat[, pvalcol] < pvalue_threshold & abs(counts_dat[, new_contrast_label]) > log2fc_threshold]
    # print(filtered_genes)
    repeated_column <- rep(contrast, length(filtered_genes))
    ## If param empty, fill it with default value.
    if (length(add_deg_columns) == 0) {
      add_deg_columns <- c("FC", "logFC", "tstat", "pval", "adjpval")
    } else if (all(add_deg_columns == "none")) {
      new_df <- data.frame(filtered_genes, repeated_column)
      names(new_df) <- c(feature_id_colname, "Contrast")
    } else {
      add_deg_columns <- setdiff(add_deg_columns, "none")
      out_columns <- paste(contrast, add_deg_columns, sep = "_")
      deg <- counts_dat[, c(feature_id_colname, out_columns)]
      names(deg)[1] <- feature_id_colname
      new_df <- data.frame(filtered_genes, repeated_column) %>% dplyr::left_join(deg, by = c("filtered_genes" = feature_id_colname))
      names(new_df) <- c(feature_id_colname, "Contrast", add_deg_columns)
    }

    df_out1 <- new_df
    df_outs[[contrast]] <- df_out1

    ### PH: END Build table for Volcano plot


    ### PH: START Make plot - Can we use Enhanced volcano function from other template to make figure instead of ggplot shown here

    print(paste0("Total number of genes included in volcano plot: ", nrow(counts_dat)))
    ## special nudge/repel of specific genes
    if (DisplaceGeneLabels) {
      gn <- trimws(unlist(strsplit(c(GeneListSpecialLabelDisplacement), ",")), which = c("both"))
      ind_gn <- match(gn, counts_dat$Gene[gene_list_ind]) # get the indices of the listed genes
      nudge_x_all <- rep(c(0.2), length(counts_dat$Gene[gene_list_ind]))
      nudge_y_all <- rep(c(0.2), length(counts_dat$Gene[gene_list_ind]))
      nudge_x_all[ind_gn] <- c(SpecialLabelDisplacementXAxis)
      nudge_y_all[ind_gn] <- c(SpecialLabelDisplacementYAxis)
    } else {
      nudge_x_all <- LabelXadj
      nudge_y_all <- LabelYadj
    }



    # set plot parameters
    if (UseDefaultYAxisLimit) {
      negative_log10_p_values <- -log10(counts_dat[, pvalcol])
      ymax <- ceiling(max(negative_log10_p_values[is.finite(negative_log10_p_values)]))
    } else {
      ymax <- YAxisLimit
    }
    if (UseDefaultXAxisLimit) {
      xmax1 <- ceiling(max(counts_dat[, lfccol]))
      xmax2 <- ceiling(max(-counts_dat[, lfccol]))
      xmax <- max(xmax1, xmax2)
    } else {
      xmax <- XAxisLimit
    }


    grm <- counts_dat[, c(new_contrast_label, pvalcol)]
    grm[, "neglogpval"] <- -log10(counts_dat[, pvalcol])
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
      geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), color = ColorofLogFoldChangeThresholdLine, alpha = 1.0) +
      geom_hline(yintercept = -log10(pvalue_threshold), color = ColorofPValueThresholdLine, alpha = 1.0) +
      geom_point(
        data = grm[counts_dat[, pvalcol] < pvalue_threshold, ],
        color = ColorofGenesMeetingOnlyPValueThreshold,
        size = PointSize
      ) +
      geom_point(
        data = grm[counts_dat[, pvalcol] < pvalue_threshold & abs(grm[, "FC"]) > log2fc_threshold, ],
        color = colorforgenesmeetingpvalueandfoldchangethresholds,
        size = PointSize
      ) +
      geom_text_repel(
        data = grm[gene_list_ind, ],
        label = counts_dat$Gene[gene_list_ind],
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
