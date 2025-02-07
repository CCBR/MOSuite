#' Make a heatmap
#'
#' @inheritParams filter_counts
#'
#' @return heatmap ggproto object
#' @keywords internal
#'
plot_corr_heatmap <- function(counts_dat,
                              sample_metadata,
                              sample_id_colname = NULL,
                              feature_id_colname = NULL,
                              group_colname = "Group",
                              label_colname = "Label",
                              color_values = c(
                                "#5954d6", "#e1562c", "#b80058", "#00c6f8", "#d163e6", "#00a76c",
                                "#ff9287", "#008cf9", "#006e00", "#796880", "#FFA500", "#878500"
                              )) {
  abort_packages_not_installed("amap", "ComplexHeatmap", "dendsort")

  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }

  ## Annotate
  rownames(sample_metadata) <- sample_metadata[[label_colname]]
  annoVal <- lapply(group_colname, function(x) {
    # TODO this only works on dataframes, not tibbles
    out <- as.factor(sample_metadata[, x]) %>% levels()
    # names(out)=x
    return(out)
  }) %>% unlist()
  col <- color_values[1:length(annoVal)]
  names(col) <- annoVal

  cols <- lapply(group_colname, function(x) {
    ax <- as.factor(sample_metadata[, x]) %>% levels()
    out <- col[ax]
    return(out)
  })
  names(cols) <- (group_colname)

  anno <- ComplexHeatmap::columnAnnotation(
    df = sample_metadata[, group_colname, drop = F],
    col = cols
  )


  ## Create Correlation Matrix

  old <- sample_metadata[[sample_id_colname]]
  new <- sample_metadata[[label_colname]]
  names(old) <- new
  counts_dat %<>% dplyr::rename(tidyselect::any_of(old))
  if (!is.null(feature_id_colname) && feature_id_colname %in% colnames(counts_dat)) {
    counts_dat %<>%
      tibble::column_to_rownames(var = feature_id_colname)
  }

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
  lgd <- ComplexHeatmap::Legend(new.palette(20), title = "Correlation", title_position = "lefttop-rot")
  hm <- ComplexHeatmap::Heatmap(m,
    heatmap_legend_param = list(
      title = "Correlation",
      title_position = "leftcenter-rot"
    ),
    cluster_rows = dend,
    cluster_columns = dend,
    top_annotation = anno,
    row_names_gp = grid::gpar(fontsize = 15),
    column_names_gp = grid::gpar(fontsize = 15),
    col = new.palette(20)
  )

  return(hm)
}
