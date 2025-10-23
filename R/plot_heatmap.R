#' Plot correlation heatmap
#'
#' @inherit moo_counts description
#'
#'
#' @param moo_counts counts dataframe or `multiOmicDataSet` containing `count_type` & `sub_count_type` in the counts
#'   slot
#' @param ... arguments forwarded to method [plot_corr_heatmap_dat]
#'
#' @export
#' @returns heatmap from `ComplexHeatmap::Heatmap()`
#' @examples
#' # plot correlation heatmap for a counts slot in a multiOmicDataset Object
#' moo <- multiOmicDataSet(
#'   sample_metadata = as.data.frame(nidap_sample_metadata),
#'   anno_dat = data.frame(),
#'   counts_lst = list("raw" = as.data.frame(nidap_raw_counts))
#' )
#' p <- plot_corr_heatmap(moo, count_type = "raw")
#'
#' # plot correlation heatmap for a counts dataframe
#' plot_corr_heatmap(
#'   moo@counts$raw,
#'   sample_metadata = moo@sample_meta,
#'   sample_id_colname = "Sample",
#'   feature_id_colname = "Gene",
#'   group_colname = "Group",
#'   label_colname = "Label"
#' )
#' @details
#'
#' # Methods
#'
#' | link to docs  | class  |
#' |---|---|
#' | [plot_corr_heatmap_moo] | `multiOmicDataSet` |
#' | [plot_corr_heatmap_dat] | `data.frame`       |
#'
#' ## Method Usage
#'
#' ```
#' # multiOmicDataSet
#' plot_corr_heatmap(moo_counts,
#'   count_type,
#'   sub_count_type = NULL,
#'   ...)
#'
#' # dataframe
#' plot_corr_heatmap(moo_counts,
#'   sample_metadata,
#'   sample_id_colname = NULL,
#'   feature_id_colname = NULL,
#'   group_colname = "Group",
#'   label_colname = "Label",
#'   color_values = c(
#'     "#5954d6", "#e1562c", "#b80058", "#00c6f8", "#d163e6", "#00a76c",
#'     "#ff9287", "#008cf9", "#006e00", "#796880", "#FFA500", "#878500"
#'   ))
#' ```
#'
#' @rdname plot_corr_heatmap
#' @family plotters
#' @family heatmaps
#' @keywords plotters
#' @family moo methods
plot_corr_heatmap <- S7::new_generic("plot_corr_heatmap", "moo_counts")

#' Plot correlation heatmap for multiOmicDataSet
#'
#' @param moo_counts `multiOmicDataSet` containing `count_type` & `sub_count_type` in the counts slot
#' @param count_type the type of counts to use. Must be a name in the counts slot (`names(moo@counts)`).
#' @param sub_count_type used if `count_type` is a list in the counts slot: specify the sub count type within the list.
#'   Must be a name in `names(moo@counts[[count_type]])`.
#' @param ... arguments forwarded to method [plot_corr_heatmap_dat]
#'
#' @name plot_corr_heatmap_moo
#' @method plot_corr_heatmap multiOmicDataSet
#' @seealso [plot_corr_heatmap] generic
#' @family plotters for multiOmicDataSets
S7::method(plot_corr_heatmap, multiOmicDataSet) <- function(moo_counts,
                                                            count_type,
                                                            sub_count_type = NULL,
                                                            ...) {
  counts_dat <- extract_counts(moo_counts, count_type, sub_count_type)
  return(plot_corr_heatmap(counts_dat, sample_metadata = moo_counts@sample_meta, ...))
}

#' Plot correlation heatmap for counts dataframe
#'
#' @param moo_counts counts dataframe (**Required**)
#' @param sample_metadata sample metadata as a data frame or tibble (**Required**)
#' @param sample_id_colname The column from the sample metadata containing the sample names. The names in this column
#'   must exactly match the names used as the sample column names of your input Counts Matrix. (Default: `NULL` - first
#'   column in the sample metadata will be used.)
#' @param feature_id_colname The column from the counts data containing the Feature IDs (Usually Gene or Protein ID).
#'   This is usually the first column of your input Counts Matrix. Only columns of Text type from your input Counts
#'   Matrix will be available to select for this parameter. (Default: `NULL` - first column in the counts matrix will be
#'   used.)
#' @param group_colname The column from the sample metadata containing the sample group information. This is usually a
#'   column showing to which experimental treatments each sample belongs (e.g. WildType, Knockout, Tumor, Normal,
#'   Before, After, etc.).
#' @param label_colname The column from the sample metadata containing the sample labels as you wish them to appear in
#'   the plots produced by this template. This can be the same Sample Names Column. However, you may desire different
#'   labels to display on your figure (e.g. shorter labels are sometimes preferred on plots). In that case, select the
#'   column with your preferred Labels here. The selected column should contain unique names for each sample. (Default:
#'   `NULL` -- `sample_id_colname` will be used.)
#' @param color_values vector of colors as hex values or names recognized by R
#'
#' @name plot_corr_heatmap_dat
#' @seealso [plot_corr_heatmap] generic
#' @family plotters for counts dataframes
S7::method(plot_corr_heatmap, S7::class_data.frame) <- function(moo_counts,
                                                                sample_metadata,
                                                                sample_id_colname = NULL,
                                                                feature_id_colname = NULL,
                                                                group_colname = "Group",
                                                                label_colname = "Label",
                                                                color_values = c(
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
                                                                )) {
  abort_packages_not_installed("amap", "ComplexHeatmap", "dendsort")
  counts_dat <- moo_counts
  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (!is.null(feature_id_colname) && feature_id_colname %in% colnames(counts_dat)) {
    counts_dat <- counts_dat %>%
      tibble::column_to_rownames(var = feature_id_colname)
  }
  # drop non-numeric columns
  counts_dat <- counts_dat %>% dplyr::select(tidyselect::where(is.numeric))

  ## Annotate
  rownames(sample_metadata) <- sample_metadata[[label_colname]]
  annoVal <- lapply(group_colname, function(x) {
    # TODO this only works on dataframes, not tibbles
    out <- as.factor(sample_metadata[, x]) %>% levels()
    # names(out)=x
    return(out)
  }) %>%
    unlist()
  col <- color_values[seq_along(annoVal)]
  names(col) <- annoVal

  cols <- lapply(group_colname, function(x) {
    ax <- as.factor(sample_metadata[, x]) %>% levels()
    out <- col[ax]
    return(out)
  })
  names(cols) <- (group_colname)

  anno <- ComplexHeatmap::columnAnnotation(df = sample_metadata[, group_colname, drop = FALSE], col = cols)


  ## Create Correlation Matrix

  old <- sample_metadata[[sample_id_colname]]
  new <- sample_metadata[[label_colname]]
  names(old) <- new
  counts_dat <- counts_dat %>% dplyr::rename(tidyselect::any_of(old))

  mat <- as.matrix(counts_dat)
  tcounts <- t(mat)


  ## calculate correlation
  d <- amap::Dist(tcounts, method = "correlation", diag = TRUE)
  m <- as.matrix(d)

  ## create dendogram
  dend <- d %>%
    stats::hclust(method = "average") %>%
    stats::as.dendrogram() %>%
    dendsort::dendsort() %>%
    rev()

  ### plot
  new.palette <- grDevices::colorRampPalette(c("blue", "green", "yellow"))
  # lgd <- ComplexHeatmap::Legend(new.palette(20),
  #                               title = "Correlation",
  #                               title_position = "lefttop-rot")
  hm <- ComplexHeatmap::Heatmap(
    m,
    heatmap_legend_param = list(title = "Correlation", title_position = "leftcenter-rot"),
    cluster_rows = dend,
    cluster_columns = dend,
    top_annotation = anno,
    row_names_gp = grid::gpar(fontsize = 15),
    column_names_gp = grid::gpar(fontsize = 15),
    col = new.palette(20)
  )

  return(hm)
}

#' Plot expression heatmap
#'
#' @inherit moo_counts description
#'
#' @param moo_counts counts dataframe or `multiOmicDataSet` containing `count_type` & `sub_count_type` in the counts
#'   slot
#' @param ... arguments forwarded to method
#'
#' @export
#' @returns heatmap from `ComplexHeatmap::Heatmap()`
#'
#' @examples
#' # plot expression heatmap for a counts slot in a multiOmicDataset Object
#' moo <- multiOmicDataSet(
#'   sample_metadata = as.data.frame(nidap_sample_metadata),
#'   anno_dat = data.frame(),
#'   counts_lst = list(
#'     "raw" = nidap_raw_counts,
#'     "norm" = list(
#'       "voom" = as.data.frame(nidap_norm_counts)
#'     )
#'   )
#' )
#' p <- plot_expr_heatmap(moo, count_type = "norm", sub_count_type = "voom")
#'
#' # customize the plot
#' plot_expr_heatmap(moo,
#'   count_type = "norm", sub_count_type = "voom",
#'   top_genes_by_variance_to_include = 100
#' )
#'
#' # plot expression heatmap for a counts dataframe
#' counts_dat <- moo@counts$norm$voom
#' plot_expr_heatmap(
#'   counts_dat,
#'   sample_metadata = nidap_sample_metadata,
#'   sample_id_colname = "Sample",
#'   feature_id_colname = "Gene",
#'   group_colname = "Group",
#'   label_colname = "Label",
#'   top_genes_by_variance_to_include = 100
#' )
#'
#' @details
#'
#' The samples (i.e. the columns) are clustered in an unsupervised fashion based
#' on how similar their expression profiles are across the included genes. This
#' can help identify samples that are non clustering with their group as you
#' might expect based on the experimental design.
#'
#' By default, the top 500 genes by variance are used, as these are
#' generally going to include those genes that most distinguish your samples
#' from one another. You can change this as well as many other parameters about
#' this heatmap if you explore the advanced options (see [plot_expr_heatmap_dat]).
#'
#' # Methods
#'
#' | link to docs  | class  |
#' |---|---|
#' | [plot_expr_heatmap_moo] | `multiOmicDataSet` |
#' | [plot_expr_heatmap_dat] | `data.frame`       |
#'
#' ## Method Usage
#'
#' ```
#' # multiOmicDataset
#' plot_expr_heatmap(moo_counts,
#'   count_type,
#'   sub_count_type = NULL,
#'   ...)
#'
#' # dataframe
#' plot_expr_heatmap(moo_counts,
#'   sample_metadata,
#'   sample_id_colname = NULL,
#'   feature_id_colname = NULL,
#'   group_colname = "Group",
#'   label_colname = NULL,
#'   samples_to_include = NULL,
#'   color_values = c(
#'     "#5954d6", "#e1562c", "#b80058", "#00c6f8", "#d163e6", "#00a76c",
#'     "#ff9287", "#008cf9", "#006e00", "#796880", "#FFA500", "#878500"
#'   ),
#'   include_all_genes = FALSE,
#'   filter_top_genes_by_variance = TRUE,
#'   top_genes_by_variance_to_include = 500,
#'   specific_genes_to_include_in_heatmap = "None",
#'   cluster_genes = TRUE,
#'   gene_distance_metric = "correlation",
#'   gene_clustering_method = "average",
#'   display_gene_dendrograms = TRUE,
#'   display_gene_names = FALSE,
#'   center_and_rescale_expression = TRUE,
#'   cluster_samples = FALSE,
#'   arrange_sample_columns = TRUE,
#'   order_by_gene_expression = FALSE,
#'   gene_to_order_columns = " ",
#'   gene_expression_order = "low_to_high",
#'   smpl_distance_metric = "correlation",
#'   smpl_clustering_method = "average",
#'   display_smpl_dendrograms = TRUE,
#'   reorder_dendrogram = FALSE,
#'   reorder_dendrogram_order = c(),
#'   display_sample_names = TRUE,
#'   group_columns = c("Group", "Replicate", "Batch"),
#'   assign_group_colors = FALSE,
#'   assign_color_to_sample_groups = c(),
#'   group_colors = c("indigo", "carrot", "lipstick", "turquoise", "lavender",
#'     "jade", "coral", "azure", "green", "rum", "orange", "olive"),
#'   heatmap_color_scheme = "Default",
#'   autoscale_heatmap_color = TRUE,
#'   set_min_heatmap_color = -2,
#'   set_max_heatmap_color = 2,
#'   aspect_ratio = "Auto",
#'   legend_font_size = 10,
#'   gene_name_font_size = 4,
#'   sample_name_font_size = 8,
#'   display_numbers = FALSE)
#' ```
#'
#' @rdname plot_expr_heatmap
#' @family plotters
#' @family heatmaps
#' @keywords plotters
#' @family moo methods
plot_expr_heatmap <- S7::new_generic("plot_expr_heatmap", "moo_counts")


#' Plot expression heatmap for multiOmicDataSet
#'
#' @param moo_counts `multiOmicDataSet` containing `count_type` & `sub_count_type` in the counts slot
#' @param count_type the type of counts to use. Must be a name in the counts slot (`names(moo@counts)`).
#' @param sub_count_type used if `count_type` is a list in the counts slot: specify the sub count type within the list.
#'   Must be a name in `names(moo@counts[[count_type]])`.
#' @param ... arguments forwarded to method: [plot_expr_heatmap_dat]
#'
#' @name plot_expr_heatmap_moo
#' @method plot_expr_heatmap multiOmicDataSet
#' @seealso [plot_expr_heatmap] generic
#' @family plotters for multiOmicDataSets
S7::method(plot_expr_heatmap, multiOmicDataSet) <- function(moo_counts,
                                                            count_type,
                                                            sub_count_type = NULL,
                                                            ...) {
  counts_dat <- extract_counts(moo_counts, count_type, sub_count_type)
  return(plot_expr_heatmap(counts_dat, sample_metadata = moo_counts@sample_meta, ...))
}


#' Plot expression heatmap for counts dataframe
#'
#' @param moo_counts counts dataframe (**Required**)
#' @param sample_metadata sample metadata as a data frame or tibble  (**Required**)
#' @param sample_id_colname The column from the sample metadata containing the sample names. The names in this column
#'   must exactly match the names used as the sample column names of your input Counts Matrix. (Default: `NULL` - first
#'   column in the sample metadata will be used.)
#' @param feature_id_colname The column from the counts dataa containing the Feature IDs (Usually Gene or Protein ID).
#'   This is usually the first column of your input Counts Matrix. Only columns of Text type from your input Counts
#'   Matrix will be available to select for this parameter. (Default: `NULL` - first column in the counts matrix will be
#'   used.)
#' @param group_colname The column from the sample metadata containing the sample group information. This is usually a
#'   column showing to which experimental treatments each sample belongs (e.g. WildType, Knockout, Tumor, Normal,
#'   Before, After, etc.).
#' @param label_colname The column from the sample metadata containing the sample labels as you wish them to appear in
#'   the plots produced by this template. This can be the same Sample Names Column. However, you may desire different
#'   labels to display on your figure (e.g. shorter labels are sometimes preferred on plots). In that case, select the
#'   column with your preferred Labels here. The selected column should contain unique names for each sample. (Default:
#'   `NULL` -- `sample_id_colname` will be used.)
#' @param color_values vector of colors as hex values or names recognized by R
#' @param samples_to_include Which samples would you like to include? Usually, you will choose all sample columns, or
#'   you could choose to remove certain samples. Samples excluded here will be removed in this step and from further
#'   analysis downstream of this step. (Default: `NULL` - all sample IDs in `moo@sample_meta` will be used.)
#' @param include_all_genes Set to TRUE if all genes are to be included. Set to FALSE if you want to filter genes by
#'   variance and/or provide a list of specific genes that will appear in the heatmap.
#' @param filter_top_genes_by_variance Set to TRUE if you want to only include the top genes by variance. Set to FALSE
#'   if you do not want to filter genes by variance.
#' @param top_genes_by_variance_to_include The number of genes to include if filtering genes by variance. This parameter
#'   is ignored if "Filter top genes by variance" is set to FALSE.
#' @param specific_genes_to_include_in_heatmap Enter the gene symbols to be included in the heatmap, with each gene
#'   symbol separated with a space from the others. Alternatively, paste in a column of gene names from any spreadsheet
#'   application. This parameter is ignored if "Include all genes" is set to TRUE.
#' @param cluster_genes Choose whether to cluster the rows (genes). If TRUE, rows will have clustering applied. If
#'   FALSE, clustering will not be applied to rows.
#' @param gene_distance_metric Distance metric to be used in clustering genes. (TODO document options)
#' @param gene_clustering_method Clustering method metric to be used in clustering samples. (TODO document options)
#' @param display_gene_dendrograms Set to TRUE to show gene dendrograms. Set to FALSE to hide dendrograms.
#' @param display_gene_names Set to TRUE to display gene names on the right side of the heatmap. Set to FALSE to hide
#'   gene names.
#' @param center_and_rescale_expression Center and rescale expression for each gene across all included samples.
#' @param cluster_samples Choose whether to cluster the columns (samples). If TRUE, columns will have clustering
#'   applied. If FALSE, clustering will not be applied to columns.
#' @param arrange_sample_columns If TRUE, arranges columns by annotation groups. If FALSE, and "Cluster Samples" is
#'   FALSE, samples will appear in the order of input (samples to include)
#' @param order_by_gene_expression If TRUE, set gene name below and direction for ordering
#' @param gene_to_order_columns Gene to order columns by expression levels
#' @param gene_expression_order Choose direction for gene order
#' @param smpl_distance_metric Distance metric to be used in clustering samples.  (TODO document options)
#' @param smpl_clustering_method Clustering method to be used in clustering samples.  (TODO document options)
#' @param display_smpl_dendrograms Set to TRUE to show sample dendrograms. Set to FALSE to hide dendrogram.
#' @param reorder_dendrogram If TRUE, set the order of the dendrogram (below)
#' @param reorder_dendrogram_order Reorder the samples (columns) of the dendrogram by name, e.g.
#'   “sample2”,“sample3",“sample1".
#' @param display_sample_names Set to TRUE if you want sample names to be displayed on the plot. Set to FALSE to hide
#'   sample names.
#' @param group_columns Columns containing the sample groups for annotation tracks
#' @param assign_group_colors If TRUE, set the groups assigned colors (below)
#' @param assign_color_to_sample_groups Enter each sample to color in the format: group_name: color This parameter is
#'   ignored if "Assign Colors" is set to FALSE.
#' @param group_colors Set group annotation colors.
#' @param heatmap_color_scheme color scheme (TODO document options)
#' @param autoscale_heatmap_color Set to TRUE to autoscale the heatmap colors between the maximum and minimum heatmap
#'   color parameters. If FALSE, set the heatmap colors between "Set max heatmap color" and "Set min heatmap color"
#'   (below).
#' @param set_min_heatmap_color If Autoscale heatmap color is set to FALSE, set the minimum heatmap z-score value
#' @param set_max_heatmap_color If Autoscale heatmap color is set to FALSE, set the maximum heatmap z-score value.
#' @param aspect_ratio Set figure Aspect Ratio. Ratio refers to entire figure including legend. If set to Auto figure
#'   size is based on number of rows and columns form counts matrix. default - Auto
#' @param legend_font_size Set Font size for figure legend. Default is 10.
#' @param gene_name_font_size Font size for gene names. If you don't want gene labels to show, toggle "Display Gene
#'   Names" below to FALSE
#' @param sample_name_font_size Font size for sample names. If you don't want to display samples names, toggle "Display
#'   sample names" (below) to FALSE
#' @param display_numbers Setting to FALSE (default) will not display numerical value of heat on heatmap. Set to TRUE if
#'   you want to see these numbers on the plot.
#'
#' @name plot_expr_heatmap_dat
#' @seealso [plot_expr_heatmap] generic
#' @family heatmaps
#' @family plotters for counts dataframes
S7::method(plot_expr_heatmap, S7::class_data.frame) <- function(moo_counts,
                                                                sample_metadata,
                                                                sample_id_colname = NULL,
                                                                feature_id_colname = NULL,
                                                                group_colname = "Group",
                                                                label_colname = NULL,
                                                                samples_to_include = NULL,
                                                                color_values = c(
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
                                                                include_all_genes = FALSE,
                                                                filter_top_genes_by_variance = TRUE,
                                                                top_genes_by_variance_to_include = 500,
                                                                specific_genes_to_include_in_heatmap = "None",
                                                                cluster_genes = TRUE,
                                                                gene_distance_metric = "correlation",
                                                                gene_clustering_method = "average",
                                                                display_gene_dendrograms = TRUE,
                                                                display_gene_names = FALSE,
                                                                center_and_rescale_expression = TRUE,
                                                                cluster_samples = FALSE,
                                                                arrange_sample_columns = TRUE,
                                                                order_by_gene_expression = FALSE,
                                                                gene_to_order_columns = " ",
                                                                gene_expression_order = "low_to_high",
                                                                smpl_distance_metric = "correlation",
                                                                smpl_clustering_method = "average",
                                                                display_smpl_dendrograms = TRUE,
                                                                reorder_dendrogram = FALSE,
                                                                reorder_dendrogram_order = c(),
                                                                display_sample_names = TRUE,
                                                                group_columns = c("Group", "Replicate", "Batch"),
                                                                assign_group_colors = FALSE,
                                                                assign_color_to_sample_groups = c(),
                                                                group_colors = c(
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
                                                                ),
                                                                heatmap_color_scheme = "Default",
                                                                autoscale_heatmap_color = TRUE,
                                                                set_min_heatmap_color = -2,
                                                                set_max_heatmap_color = 2,
                                                                aspect_ratio = "Auto",
                                                                legend_font_size = 10,
                                                                gene_name_font_size = 4,
                                                                sample_name_font_size = 8,
                                                                display_numbers = FALSE) {
  ## This function uses pheatmap to draw a heatmap, scaling first by rows
  ## (with samples in columns and genes in rows)
  Gene <- NULL
  # TODO support tibbles; currently these must be dataframes
  counts_dat <- as.data.frame(moo_counts)
  sample_metadata <- as.data.frame(sample_metadata)


  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(counts_dat)[1]
  }
  if (is.null(label_colname)) {
    label_colname <- sample_id_colname
  }
  if (is.null(samples_to_include)) {
    samples_to_include <- sample_metadata %>% dplyr::pull(sample_id_colname)
  }

  ## --------------- ##
  ## Error Messages ##
  ## -------------- ##

  if (include_all_genes == TRUE && filter_top_genes_by_variance == TRUE) {
    stop("ERROR: Choose only one of 'Include all genes' or 'Filter top genes by variance' as TRUE")
  }

  if (any(
    all(cluster_samples == TRUE, arrange_sample_columns == TRUE),
    all(
      arrange_sample_columns == TRUE,
      order_by_gene_expression == TRUE
    ),
    all(arrange_sample_columns == TRUE, cluster_samples == TRUE),
    all(
      cluster_samples == FALSE,
      arrange_sample_columns == FALSE,
      order_by_gene_expression == FALSE
    )
  )) {
    stop(
      "ERROR: Choose only one of 'Cluster Samples', 'Arrange sample columns', or 'order by gene expression' as TRUE"
    )
  }

  ### PH: START palette function for heatmap scale
  ## Begin pal() color palette function∂:
  pal <- function(n,
                  h = c(237, 43),
                  c = 100,
                  l = c(70, 90),
                  power = 1,
                  fixup = TRUE,
                  gamma = NULL,
                  alpha = 1,
                  ...) {
    if (n < 1L) {
      return(character(0L))
    }
    h <- rep(h, length.out = 2L)
    c <- c[1L]
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, -1, length = n)
    rval <- colorspace::hex(
      colorspace::polarLUV(
        L = l[2L] - diff(l) * abs(rval)^power[2L],
        C = c * abs(rval)^power[1L],
        H = ifelse(rval > 0, h[1L], h[2L])
      ),
      fixup = fixup,
      ...
    )
    if (!missing(alpha)) {
      alpha <- pmax(pmin(alpha, 1), 0)
      alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)),
        width = 2L,
        upper.case = TRUE
      )
      rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
  }
  ### PH: END palette function for heatmap scale

  ### PH: START SET up heatmap function for do.call
  ## Stratagy is to use Pheatmap to create heatmap then output as Complex Heatmap to add Annotations

  ## Begin doheatmap() function:
  doheatmap <- function(dat, clus, clus2, ht, rn, cn, col, dispnum) {
    col.pal <- np[[col]]
    # if (=FALSE) {
    #   col.pal = rev(col.pal)
    # }
    # Define metrics for clustering
    drows1 <- gene_distance_metric
    dcols1 <- smpl_distance_metric
    minx <- min(dat)
    maxx <- max(dat)
    if (autoscale_heatmap_color) {
      breaks <- seq(minx, maxx, length = 100)
      legbreaks <- seq(minx, maxx, length = 5)
    } else {
      breaks <- seq(set_min_heatmap_color, set_max_heatmap_color, length = 100)
      legbreaks <- seq(set_min_heatmap_color, set_max_heatmap_color, length = 5)
    }
    breaks <- sapply(breaks, signif, 4)
    legbreaks <- sapply(legbreaks, signif, 4)
    # Run cluster method using
    hcrow <- stats::hclust(stats::dist(dat), method = gene_clustering_method)
    # hc <- stats::hclust(stats::dist(t(dat)), method = smpl_clustering_method)

    if (FALSE) {
      sort_hclust <- function(...) {
        return(stats::as.hclust(rev(dendsort::dendsort(stats::as.dendrogram(...)))))
      }
    } else {
      sort_hclust <- function(...) {
        return(stats::as.hclust(dendsort::dendsort(stats::as.dendrogram(...))))
      }
    }
    # if (clus) {
    #   colclus <- sort_hclust(hc)
    # } else {
    #   colclus <- FALSE
    # }
    if (clus2) {
      rowclus <- sort_hclust(hcrow)
    } else {
      rowclus <- FALSE
    }
    if (display_smpl_dendrograms) {
      smpl_treeheight <- 25
    } else {
      smpl_treeheight <- 0
    }
    if (display_gene_dendrograms) {
      gene_treeheight <- 25
    } else {
      gene_treeheight <- 0
    }
    hm.parameters <- list(
      dat,
      color = col.pal,
      legend_breaks = legbreaks,
      legend = TRUE,
      scale = "none",
      treeheight_col = smpl_treeheight,
      treeheight_row = gene_treeheight,
      kmeans_k = NA,
      breaks = breaks,
      display_numbers = dispnum,
      number_color = "black",
      fontsize_number = 8,
      height = 80,
      cellwidth = NA,
      cellheight = NA,
      fontsize = legend_font_size,
      fontsize_row = gene_name_font_size,
      fontsize_col = sample_name_font_size,
      show_rownames = rn,
      show_colnames = cn,
      cluster_rows = rowclus,
      cluster_cols = clus,
      clustering_distance_rows = drows1,
      clustering_distance_cols = dcols1,
      annotation_col = annotation_col,
      annotation_colors = annot_col,
      labels_col = labels_col
    )
    # mat <- t(dat)
    callback <- function(hc, mat) {
      dend <- rev(dendsort::dendsort(stats::as.dendrogram(hc)))
      if (reorder_dendrogram == TRUE) {
        dend <- dend %>% dendextend::rotate(reorder_dendrogram_order)
      } else {
        dend <- dend %>% dendextend::rotate(c(1:stats::nobs(dend)))
      }
      return(stats::as.hclust(dend))
    }
    ### PH: END SET up heatmap function for do.call

    ## Make Heatmap
    return(do.call(ComplexHeatmap::pheatmap, c(hm.parameters, list(clustering_callback = callback))))
  }
  # End doheatmap() function.

  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##

  ### PH: START  Build different color spectra options for heatmap:
  np0 <- pal(100)
  np1 <- colorspace::diverge_hcl(100,
    c = 100,
    l = c(30, 80),
    power = 1
  ) # Blue to Red
  np2 <- colorspace::heat_hcl(
    100,
    c = c(80, 30),
    l = c(30, 90),
    power = c(1 / 5, 2)
  ) # Red to Vanilla
  np3 <- rev(colorspace::heat_hcl(
    100,
    h = c(0, -100),
    c = c(40, 80),
    l = c(75, 40),
    power = 1
  )) # Violet to Pink
  np4 <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(100)) # Red to yellow to blue
  np5 <- grDevices::colorRampPalette(c("steelblue", "white", "red"))(100) # Steelblue to White to Red

  ## Gather list of color spectra and give them names for the GUI to show.
  np <- list(np0, np1, np2, np3, np4, np5)
  names(np) <- c(
    "Default",
    "Blue to Red",
    "Red to Vanilla",
    "Violet to Pink",
    "Bu Yl Rd",
    "Bu Wt Rd"
  )

  ### PH: END  Build different color spectra options for heatmap:


  ### PH: START  Build Counts Table for HM

  ##############
  ### Select Samples
  ##############

  ## Parse input counts matrix. Subset by samples.
  df1 <- counts_dat
  # Swap out Gene Name column name, if it's not 'Gene'.
  # TODO: refactor to avoid renaming gene column
  if (feature_id_colname != "Gene") {
    # Drop original Gene column
    df1 <- df1[, !(colnames(df1) %in% c("Gene"))]
    # Rename column to Gene
    colnames(df1)[which(colnames(df1) == feature_id_colname)] <- "Gene"
  }
  # Build new counts matrix containing only sample subset chosen by user.
  df.orig <- df1
  df <- df.orig %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise_all(mean)
  df.mat <- df[, (colnames(df) != "Gene")] %>% as.data.frame()
  # df %>% dplyr::mutate(Gene = stringr::str_replace_all(Gene, "_", " ")) -> df
  row.names(df.mat) <- df$Gene
  rownames(df.mat) <- stringr::str_wrap(rownames(df.mat), 30) # for really long geneset names
  df.mat <- as.data.frame(df.mat)

  ##############
  ## Subset counts matrix by genes.
  ##############

  # Toggle to include all genes in counts matrix (in addition to any user-submitted gene list).
  if (include_all_genes == FALSE) {
    # Add user-submitted gene list (optional).
    genes_to_include_parsed <- c()
    genes_to_include_parsed <- strsplit(specific_genes_to_include_in_heatmap, " ")[[1]]
    # genes_to_include_parsed = gsub("_"," ",genes_to_include_parsed)
    df.final.extra.genes <- df.mat[genes_to_include_parsed, ]

    # filter all genes by variance + user-submitted gene list
    if (filter_top_genes_by_variance == TRUE) {
      df.final <- as.matrix(df.mat)
      var <- matrixStats::rowVars(df.final)
      df <- as.data.frame(df.final)
      rownames(df) <- rownames(df.final)
      df.final <- df
      df.final$var <- var
      df.final <- df.final %>%
        tibble::rownames_to_column("Gene")
      df.final <- df.final %>%
        dplyr::arrange(dplyr::desc(var))
      df.final.extra.genes <- dplyr::filter(df.final, Gene %in% genes_to_include_parsed)
      df.final <- df.final[1:top_genes_by_variance_to_include, ]
      df.final <- df.final[stats::complete.cases(df.final), ]
      # Rbind user gene list to variance-filtered gene list and deduplicate.
      df.final <- rbind(df.final, df.final.extra.genes)
      df.final <- df.final[!duplicated(df.final), ]
      rownames(df.final) <- df.final$Gene
      df.final$Gene <- NULL
      df.final$var <- NULL
    } else {
      # filter ONLY user-provided gene list.
      df.final <- df.final.extra.genes
      df.final <- df.final[!duplicated(df.final), ]
      # Order genes in heatmap by user-submitted order of gene names.
      df.final <- df.final[genes_to_include_parsed, ]
      # df.final$Gene <- NULL
    }
  } else {
    df.final <- df.mat
    df.final$Gene <- NULL
  }

  ##############
  ## Center and Rescale Counts
  ##############
  ## Optionally apply centering and rescaling (default TRUE).
  if (center_and_rescale_expression == TRUE) {
    tmean.scale <- t(scale(t(df.final)))
    tmean.scale <- tmean.scale[!is.infinite(rowSums(tmean.scale)), ]
    tmean.scale <- stats::na.omit(tmean.scale)
  } else {
    tmean.scale <- df.final
  }

  ##############
  ## Order rows by Gene Expression
  ##############
  if (order_by_gene_expression == TRUE) {
    gene_to_order_columns <- gsub(" ", "", gene_to_order_columns)
    if (gene_expression_order == "low_to_high") {
      tmean.scale <- tmean.scale[, order(tmean.scale[gene_to_order_columns, ])] # order from low to high
    } else {
      tmean.scale <- tmean.scale[, order(-tmean.scale[gene_to_order_columns, ])] # order from high to low
    }
  }

  df.final <- as.data.frame(tmean.scale)
  ### PH: END  Build Counts Table for HM


  ### PH: START  Build Annotation Columns

  ## Parse input sample metadata and add annotation tracks to top of heatmap.
  annot <- sample_metadata
  # Filter to only samples user requests.
  annot <- annot %>% dplyr::filter(.data[[sample_id_colname]] %in% samples_to_include)

  # Arrange sample options.
  if (arrange_sample_columns) {
    annot <- annot[match(samples_to_include, annot[[sample_id_colname]]), ]
    for (x in group_columns) {
      annot[, x] <- factor(annot[, x], levels = unique(annot[, x]))
    }
    annot <- annot %>% dplyr::arrange_(.dots = group_columns, .by_group = TRUE)
    df.final <- df.final[, match(annot[[sample_id_colname]], colnames(df.final))]
  }


  # Build subsetted sample metadata table to use for figure.
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
  group_colors <- colorlist[group_colors]

  annotation_col <- annot %>% dplyr::select(tidyselect::all_of(group_columns))
  annotation_col <- as.data.frame(unclass(annotation_col))
  annotation_col[] <- lapply(annotation_col, factor)
  x <- length(unlist(lapply(annotation_col, levels)))
  if (x > length(group_colors)) {
    k <- x - length(group_colors)
    more_cols <- get_random_colors(k)
    group_colors <- c(group_colors, more_cols)
  }
  rownames(annotation_col) <- annot[[label_colname]]
  annot_col <- list()
  b <- 1
  i <- 1
  while (i <= length(group_columns)) {
    cnam <- group_columns[i]
    grp <- as.factor(annotation_col[, i])
    c <- b + length(levels(grp)) - 1
    col <- group_colors[b:c]
    names(col) <- levels(grp)
    assign(cnam, col)
    annot_col <- append(annot_col, mget(cnam))
    b <- c + 1
    i <- i + 1
  }

  if (assign_group_colors == TRUE) {
    colassign <- assign_color_to_sample_groups
    groupname <- c()
    groupcol <- c()
    for (i in seq_along(colassign)) {
      groupname[i] <- strsplit(colassign[i], ": ?")[[1]][1]
      groupcol[i] <- strsplit(colassign[i], ": ?")[[1]][2]
    }
    annot_col[[1]][groupname] <- groupcol
  }
  ### PH: End  Build Annotation Columns

  old <- annot[[sample_id_colname]]
  new <- annot[[label_colname]]
  names(old) <- new
  df.final <- dplyr::rename(df.final, tidyselect::any_of(old))
  labels_col <- colnames(df.final)


  ## Print number of genes to log.
  print(paste0("The total number of genes in heatmap: ", nrow(df.final)))

  ## PH: Make the heatmap.
  p <- doheatmap(
    dat = df.final,
    clus = cluster_samples,
    clus2 = cluster_genes,
    ht = 50,
    rn = display_gene_names,
    cn = display_sample_names,
    col = heatmap_color_scheme,
    dispnum = display_numbers
  )
  p@matrix_color_mapping@name <- " "
  p@matrix_legend_param$at <- as.numeric(formatC(p@matrix_legend_param$at, 2))
  p@column_title_param$gp$fontsize <- 10
  # print(p)

  ## PH: Output heatmap counts table.
  ## If user sets toggle to TRUE, return Z-scores.
  ## Else return input counts matrix by default (toggle FALSE).
  ## Returned matrix includes only genes & samples used in heatmap.
  # if (return_z_scores) {
  #   df.new <- data.frame(tmean.scale) # Convert to Z-scores.
  #   df.new %>% tibble::rownames_to_column("Gene") -> df.new
  #   return(df.new)
  # } else {
  #   df.final %>% tibble::rownames_to_column("Gene") -> df.new
  #   return(df.new)
  # }
  return(p)
}
