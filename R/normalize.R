#' Normalize counts
#'
#' @inheritParams filter_counts
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
#'   normalize(
#'     gene_names_column = "Gene",
#'     columns_to_include = c("Gene", "A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"),
#'     sample_names_column = "Sample",
#'     group_column = "Group",
#'     label_column = "Label"
#'   )
#' head(moo@counts[["norm"]][["voom"]])
normalize <- function(moo,
                      count_type = "filt",
                      gene_names_column = "Gene",
                      columns_to_include = c("Gene", "A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"),
                      sample_names_column = "Sample",
                      group_column = "Group",
                      label_column = "Label",
                      input_in_log_counts = FALSE,
                      normalization_method = "quantile",
                      samples_to_rename_manually_on_pca = c(""),
                      add_labels_to_pca = TRUE,
                      principal_component_on_x_axis_for_pca = 1,
                      principal_component_on_y_axis_for_pca = 2,
                      legend_position_for_pca = "top",
                      label_offset_x_for_pca = 2,
                      label_offset_y_for_pca = 2,
                      label_font_size_for_pca = 3,
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
  library(limma)
  library(tidyverse)
  library(edgeR)
  library(ggplot2)
  library(plotly)
  library(dplyr)
  library(RColorBrewer)
  library(colorspace)
  library(stringr)
  library(RCurl)
  library(reshape2)
  library(gridExtra)
  library(amap)
  library(lattice)
  library(gplots)
  library(gridGraphics)
  library(dendsort)
  library(ComplexHeatmap)
  library(ggrepel)

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
    x <- DGEList(counts = 2^df.filt, genes = gene_names)
  } else {
    x <- DGEList(counts = df.filt, genes = gene_names)
  }

  v <- voom(x, normalize = normalization_method)
  rownames(v$E) <- v$genes$GeneID
  as.data.frame(v$E) %>% rownames_to_column(gene_names_column) -> df.voom
  print(paste0("Total number of features included: ", nrow(df.voom)))

  ### PH: END Limma Normalization



  ### PH: START PCA Plot
  ########################
  ## PCA Plot:
  ########################
  ## same function as used from Filter Template except input is Normalized Counts table instead of Raw Data.
  ## major difference is that rawData had to be log transformed and This data does not so we will either
  ## need a parameter that instructs function to log transform data or some logic to detect log transformed data


  edf <- v$E
  tedf <- t(edf)
  tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
  tedf <- tedf[, apply(tedf, 2, var) != 0]
  pca <- prcomp(tedf, scale. = T)

  pcx <- paste0("PC", principal_component_on_x_axis_for_pca)
  pcy <- paste0("PC", principal_component_on_y_axis_for_pca)
  pca.df <- as.data.frame(pca$x) %>% dplyr::select(.data[[pcx]], .data[[pcy]])
  pca.df$group <- sample_metadata[[group_column]]
  pca.df$sample <- sample_metadata[[label_column]]
  perc.var <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  perc.var <- formatC(perc.var, format = "g", digits = 4)
  pc.x.lab <- paste0(pcx, " ", perc.var[principal_component_on_x_axis_for_pca], "%")
  pc.y.lab <- paste0(pcy, " ", perc.var[principal_component_on_y_axis_for_pca], "%")
  labelpos <- pca.df
  labelpos$mean_y <- pca.df[[pcy]] + label_offset_y_for_pca
  labelpos$mean_x <- pca.df[[pcx]] + label_offset_x_for_pca
  pca.df$xdata <- pca.df[[pcx]]
  pca.df$ydata <- pca.df[[pcy]]

  # Manual changes to sample names
  replacements <- samples_to_rename_manually_on_pca

  if (!is.null(replacements)) {
    if (replacements != c("")) {
      for (x in replacements) {
        old <- strsplit(x, ": ?")[[1]][1]
        new <- strsplit(x, ": ?")[[1]][2]
        pca.df$sample <- ifelse(pca.df$sample == old, new, pca.df$sample)
      }
    }
  }


  colorval <- colorlist[colors_for_plots]
  colorval <- unname(colorval) # remove names which affect ggplot

  if (length(unique(sample_metadata[[group_column]])) > length(colorval)) {
    ## Original color-picking code.
    k <- length(unique(sample_metadata[[group_column]])) - length(colorval)
    more_cols <- get_random_colors(k)
    colorval <- c(colorval, more_cols)
  }
  #### plot PCA
  pcaPlot <- ggplot(pca.df, aes(x = xdata, y = ydata, text = sample)) +
    geom_point(aes(color = group), text = sample, size = point_size_for_pca) +
    theme_bw() +
    theme(
      legend.position = legend_position_for_pca,
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 20),
      # panel.grid.major = element_line(size = 1),
      # axis.line=element_line(size=1),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 1
      ),
      axis.ticks = element_line(size = 1),
      legend.text = element_text(size = 18)
    ) +
    coord_fixed() +
    scale_colour_manual(values = colorval) +
    xlab(pc.x.lab) +
    ylab(pc.y.lab)

  if (add_labels_to_pca == TRUE) {
    pcaPlot <- pcaPlot +
      geom_text_repel(
        aes(label = sample, color = group),
        size = 7,
        show.legend = F,
        direction = c("both"),
        box.padding = 1.25
      )
  }

  par(mfrow = c(2, 1))
  ### PH: END PCA Plot


  ### PH: START Histogram plot
  ########################
  ## Start Histogram Plot:
  ########################
  ## same function as used from Filter Template only change I think we need is to recognize the type of Counts input to label X axis ( is this normalized or Raw Counts)
  ## my thought would be just make it a parameter to set  X axis label and we will set the defaults in the Workflow function

  df.m <- melt(edf, id.vars = c(gene_names_column))
  df.m <- dplyr::rename(df.m, sample = Var2)

  if (set_min_max_for_x_axis_for_histogram == TRUE) {
    xmin <- minimum_for_x_axis_for_histogram
    xmax <- maximum_for_x_axis_for_histogram
  } else {
    xmin <- min(df.m$value)
    xmax <- max(df.m$value)
  }

  if (color_histogram_by_group == TRUE) {
    df.m %>% mutate(colgroup = sample_metadata[sample, group_column]) -> df.m
    df.m <- df.m[complete.cases(df.m[, "colgroup"]), ]
    df.m$colgroup <- gsub("\\s", "_", df.m$colgroup)
    df.m$colgroup <- factor(df.m$colgroup, levels = unique(df.m$colgroup))
    # print(unique(df.m$sample))
    n <- length(levels(df.m$colgroup))
    cols <- colorval[1:n]
    # cols<- getourrandomcolors(n)

    histPlot <- ggplot(df.m, aes(x = value, group = sample)) +
      geom_density(aes(colour = colgroup), size = 1)
  } else {
    df.m$sample <- sample_metadata[df.m$sample, label_column]
    n <- length(unique(df.m$sample))
    cols <- getourrandomcolors(n)
    # cols=colorval[1:n]

    histPlot <- ggplot(df.m, aes(x = value, group = sample)) +
      geom_density(aes(colour = sample), size = 1)
  }
  histPlot <- histPlot +
    xlab("Normalized Counts") + ylab("Density") +
    theme_bw() +
    theme(
      legend.position = legend_position_for_histogram,
      legend.text = element_text(size = legend_font_size_for_histogram),
      legend.title = element_blank(),
      panel.background = element_blank(),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 20),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 0
      ),
      axis.line = element_line(size = .5),
      axis.ticks = element_line(size = 1)
    ) +
    # ggtitle("Frequency Histogram") +
    xlim(xmin, xmax) +
    # scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),n)) +
    scale_colour_manual(values = cols) +
    guides(linetype = guide_legend(ncol = number_of_histogram_legend_columns))

  ########################
  ### Output Figures
  ########################

  ### PH START:  Make Plots interactive.
  ## We can add logic to PCA and Histrogram to make them interactive using Plotly
  if (plot_correlation_matrix_heatmap == TRUE) {
    if (make_plots_interactive == TRUE) {
      pcaPlot1 <- (pcaPlot) %>% ggplotly(tooltip = c("sample", "group"))
      histPlot2 <- (histPlot + theme(legend.position = "none")) %>% ggplotly(tooltip = c("sample"))
      ### PH END:  Make Plots interactive.
      ## Rest of Code is NIDAP specific


      grid.newpage()
      print(pcaPlot1)
      grid.newpage()
      print(histPlot2)
    } else {
      require(gridExtra)
      corHM <- plot_heatmap(
        counts_matrix = df.filt,
        sample_metadata = sample_metadata,
        anno_col = colorval,
        anno_column = group_column,
        label_column = label_column,
        sample_names_column = sample_names_column
      )


      grid.newpage()
      print(pcaPlot)
      grid.newpage()
      print(corHM)
      grid.newpage()
      print(histPlot)
    }
  } else {
    if (make_plots_interactive == TRUE) {
      pcaPlot1 <- (pcaPlot) %>% ggplotly(tooltip = c("sample", "group"))
      histPlot2 <- (histPlot + theme(legend.position = "none")) %>% ggplotly(tooltip = "sample")

      grid.newpage()
      print(pcaPlot1)
      grid.newpage()
      print(histPlot2)
    } else {
      grid.newpage()
      print(pcaPlot)
      grid.newpage()
      print(histPlot)
    }
  }

  print("Sample columns")
  print(colnames(df.voom)[!colnames(df.voom) %in% gene_names_column])
  print("Feature Columns")
  print(colnames(anno_tbl))

  df.voom <- merge(anno_tbl, df.voom, by = gene_names_column, all.y = T)
  df.voom[, gene_names_column] <- gsub("_[0-9]+$", "", df.voom[, gene_names_column])

  if (isFALSE("norm" %in% names(moo@counts))) {
    moo@counts[["norm"]] <- list()
  }
  moo@counts[["norm"]][["voom"]] <- df.voom
  return(moo)
}
