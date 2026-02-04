#' Differential expression analysis
#'
#' @inheritParams filter_counts
#' @inheritParams batch_correct_counts
#' @inheritParams normalize_counts
#' @inheritParams option_params
#'
#' @param sub_count_type if `count_type` is a list, specify the sub count type within the list. (Default: `NULL`)
#' @param covariates_colnames The column name(s) from the sample metadata containing variable(s) of interest, such as
#'   phenotype. Most commonly this will be the same column selected for your Groups Column. Some experimental designs
#'   may require that you add additional covariate columns here.
#' @param contrast_colname The column in the metadata that contains the group variables you wish to find differential
#'   expression between. Up to 2 columns (2-factor analysis) can be used.
#' @param contrasts Specify each contrast in the format group1-group2, e.g. treated-control
#' @param covariates_colnames Columns to be used as covariates in linear modeling. Must include column from "Contrast
#'   Variable". Most commonly your covariate will be group and batch (if you have different batches in your data).
#' @param return_mean_and_sd if TRUE, return Mean and Standard Deviation of groups in addition to DEG estimates for
#'   contrast(s)
#'
#' @returns `multiOmicDataSet` with `diff` added to the `analyses` slot (i.e. `moo@analyses$diff`)
#' @export
#'
#' @family moo methods
#'
#' @examples
#' moo <- multiOmicDataSet(
#'   sample_metadata = as.data.frame(nidap_sample_metadata),
#'   anno_dat = data.frame(),
#'   counts_lst = list(
#'     "raw" = as.data.frame(nidap_raw_counts),
#'     "clean" = as.data.frame(nidap_clean_raw_counts),
#'     "filt" = as.data.frame(nidap_filtered_counts)
#'   )
#' ) |>
#'   diff_counts(
#'     count_type = "filt",
#'     sub_count_type = NULL,
#'     sample_id_colname = "Sample",
#'     feature_id_colname = "Gene",
#'     covariates_colnames = c("Group", "Batch"),
#'     contrast_colname = c("Group"),
#'     contrasts = c("B-A", "C-A", "B-C"),
#'     voom_normalization_method = "quantile",
#'   )
#' head(moo@analyses$diff)
diff_counts <- function(
  moo,
  count_type = "filt",
  sub_count_type = NULL,
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  samples_to_include = NULL,
  covariates_colnames = NULL,
  contrast_colname = NULL,
  contrasts = NULL,
  input_in_log_counts = FALSE,
  return_mean_and_sd = FALSE,
  # return_normalized_counts = FALSE,
  voom_normalization_method = "quantile",
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = "diff"
) {
  final_res <- group <- . <- NULL
  return_normalized_counts <- FALSE
  sample_metadata <- moo@sample_meta
  message(glue::glue("* differential counts"))
  # select correct counts matrix
  if (!(count_type %in% names(moo@counts))) {
    stop(glue::glue("count_type {count_type} not in moo@counts"))
  }
  counts_dat <- moo@counts[[count_type]]
  if (!is.null(sub_count_type)) {
    if (!(inherits(counts_dat, "list"))) {
      stop(
        glue::glue(
          "{count_type} counts is not a named list. To use {count_type} counts, set sub_count_type to NULL"
        )
      )
    } else if (!(sub_count_type %in% names(counts_dat))) {
      stop(
        glue::glue(
          "sub_count_type {sub_count_type} is not in moo@counts[[{count_type}]]"
        )
      )
    }
    counts_dat <- moo@counts[[count_type]][[sub_count_type]]
  }
  if (is.null(covariates_colnames)) {
    stop("covariates_colnames vector cannot be NULL")
  }
  if (is.null(contrast_colname)) {
    stop("contrast_colname cannot be NULL")
  }
  if (is.null(contrasts)) {
    stop("contrasts vector cannot be NULL")
  }
  # ensure these are vectors, not lists. needed when using cli with JSON for args
  covariates_colnames <- covariates_colnames |> unlist()
  contrast_colname <- contrast_colname |> unlist()
  contrasts <- contrasts |> unlist()

  # TODO support tibbles
  counts_dat <- counts_dat |> as.data.frame()

  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(counts_dat)[1]
  }
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##

  ### PH: START Check Rownames - from Filtering + Normalization Template
  ## create unique rownames to correctly add back Annocolumns at end of template
  counts_dat[, feature_id_colname] <- paste0(
    counts_dat[, feature_id_colname],
    "_",
    seq_len(nrow(counts_dat))
  )

  df.m <- counts_dat
  gene_names <- NULL
  gene_names$GeneID <- counts_dat[, feature_id_colname]
  ### PH: END Check Rownames - from Filtering + Normalization Template

  ### PH: START Input Data Validation - from Filtering + Normalization Template
  ### This code block does input data validation
  # Remove samples that are not in the contrast groups:
  groups <- unique(unlist(strsplit(contrasts, "-")))
  sample_metadata <- sample_metadata |>
    dplyr::filter(.data[[contrast_colname]] %in% groups)
  df.m <- df.m |>
    dplyr::select(tidyr::all_of(c(
      feature_id_colname,
      sample_metadata |> dplyr::pull(sample_id_colname)
    )))
  ### PH: END Input Data Validation - from Filtering + Normalization Template

  ####################################
  ### Computational Functions
  ################################
  ### PH: START Create Design Formula/Table
  # Put covariates in order
  covariates_colnames <- covariates_colnames[order(
    covariates_colnames != contrast_colname
  )]

  # TODO: refactor - function to sub spaces with underscores
  for (ocv in covariates_colnames) {
    sample_metadata[[ocv]] <- gsub(
      " ",
      "_",
      sample_metadata |> dplyr::pull(ocv)
    )
  }
  contrasts <- gsub(" ", "_", contrasts)
  cov <- covariates_colnames[!covariates_colnames %in% contrast_colname]

  # Combine columns if 2-factor analysis
  if (length(contrast_colname) > 1) {
    sample_metadata <- sample_metadata |>
      dplyr::mutate(
        contmerge = paste0(
          .data[[contrast_colname[1]]],
          ".",
          .data[[contrast_colname[2]]]
        )
      )
  } else {
    sample_metadata <- sample_metadata |>
      dplyr::mutate(contmerge = .data[[contrast_colname]])
  }

  contrast_var <- factor(sample_metadata$contmerge)

  ## create Design table
  if (length(cov) > 0) {
    dm.formula <- stats::as.formula(paste(
      "~0 +",
      paste(
        "contmerge",
        paste(cov, sep = "+", collapse = "+"),
        sep = "+"
      )
    ))
    design <- stats::model.matrix(dm.formula, sample_metadata)
    colnames(design) <- gsub("contmerge", "", colnames(design))
  } else {
    dm.formula <- stats::as.formula(~ 0 + contmerge)
    design <- stats::model.matrix(dm.formula, sample_metadata)
    colnames(design) <- levels(contrast_var)
  }
  ### PH: End Create Design Formula/Table

  ### PH: START Limma Normalization - Same as in Normalize Counts
  # Create DGEList object from counts - counts should not be Log scale
  if (input_in_log_counts == TRUE) {
    df_unlog <- df.m |>
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ 2^.x))
    x <- edgeR::DGEList(counts = df_unlog, genes = gene_names)
  } else {
    x <- edgeR::DGEList(counts = df.m, genes = gene_names)
  }

  # TODO add this to existing norm function & document options
  if (
    voom_normalization_method %in% c("TMM", "TMMwzp", "RLE", "upperquartile")
  ) {
    x <- edgeR::calcNormFactors(x, method = voom_normalization_method)
    rownames(x) <- x$genes$GeneID
    v <- limma::voom(x, design = design, normalize = "none")
  } else {
    v <- limma::voom(
      x,
      design = design,
      normalize = voom_normalization_method,
      save.plot = TRUE
    )
  }
  ### PH: END Limma Normalization - Same as in Normalize Counts

  ### PH: START Linear Fit and and extract df.voom table. Could be added to Limma Normalization function above with an
  ### option to run lmFit
  rownames(v$E) <- v$genes$GeneID
  # df.voom <- as.data.frame(v$E) |> tibble::rownames_to_column(feature_id_colname)
  fit <- limma::lmFit(v, design)
  cm <- limma::makeContrasts(contrasts = contrasts, levels = design)
  ### PH: END Linear Fit and and extract df.voom table.

  ### PH: START Run Contrasts (eBays) input:
  #                    -fit from LMfit
  #                    -cm from Make Contrasts
  # Run Contrasts
  fit2 <- limma::contrasts.fit(fit, cm)
  fit2 <- limma::eBayes(fit2)
  logFC <- fit2$coefficients
  colnames(logFC) <- paste(colnames(logFC), "logFC", sep = "_")
  tstat <- fit2$t
  colnames(tstat) <- paste(colnames(tstat), "tstat", sep = "_")
  FC <- 2^fit2$coefficients
  FC <- ifelse(FC < 1, -1 / FC, FC)
  colnames(FC) <- paste(colnames(FC), "FC", sep = "_")
  pvalall <- fit2$p.value
  colnames(pvalall) <- paste(colnames(pvalall), "pval", sep = "_")
  pvaladjall <- apply(pvalall, 2, function(x) {
    return(stats::p.adjust(x, "BH"))
  })
  colnames(pvaladjall) <- paste(
    colnames(fit2$coefficients),
    "adjpval",
    sep = "_"
  )

  ### PH: END Run Contrasts (eBays) input:

  ####################################
  ### Create Output Functions
  ################################
  ### PH: START Create DEG Table
  #                    -VoomObject from Limma Normalization
  #                    -pvalall from Run Contrasts (eBays)
  #                    -pvaladjall from Run Contrasts (eBays)
  #                    -FC from Run Contrasts (eBays)
  #                    -logFC from Run Contrasts (eBays)
  #                    -tstat from Run Contrasts (eBays)
  # OUTPUT: DEG Table
  if (return_mean_and_sd == TRUE) {
    tve <- t(v$E)
    mean.df <- as.data.frame(tve) |>
      tibble::rownames_to_column(sample_id_colname) |>
      dplyr::left_join(
        sample_metadata |>
          dplyr::select(tidyr::all_of(
            c(sample_id_colname, contrast_colname)
          )),
        by = sample_id_colname
      ) |>
      dplyr::rename(group = tidyr::all_of(contrast_colname)) |>
      dplyr::group_by(group) |>
      dplyr::summarise(dplyr::across(
        dplyr::where(is.numeric),
        ~ base::mean(.x)
      )) |>
      as.data.frame()
    mat_mean <- mean.df[, -c(1, 2)] |>
      as.matrix() |>
      t()
    colnames(mat_mean) <- mean.df[, 1]
    colnames(mat_mean) <- paste(colnames(mat_mean), "mean", sep = "_")
    colnames(mat_mean) <- gsub("\\.", "_", colnames(mat_mean))
    # mat_mean <- mat_mean |> as.data.frame() |> tibble::rownames_to_column(feature_id_colname)

    sd.df <- as.data.frame(tve) |>
      tibble::rownames_to_column(sample_id_colname) |>
      dplyr::left_join(
        sample_metadata |>
          dplyr::select(tidyr::all_of(
            c(sample_id_colname, contrast_colname)
          )),
        by = sample_id_colname
      ) |>
      dplyr::rename(group = tidyr::all_of(contrast_colname)) |>
      dplyr::group_by(group) |>
      dplyr::summarise(dplyr::across(
        dplyr::where(is.numeric),
        ~ stats::sd(.x)
      )) |>
      as.data.frame()
    mat_sd <- sd.df[, -c(1, 2)] |>
      as.matrix() |>
      t()
    colnames(mat_sd) <- sd.df[, 1]
    colnames(mat_sd) <- paste(colnames(mat_sd), "sd", sep = "_")
    colnames(mat_sd) <- gsub("\\.", "_", colnames(mat_sd))
    # mat_sd <- mat_sd |> as.data.frame() |> tibble::rownames_to_column(feature_id_colname)

    finalres <- purrr::map(
      list(mat_mean, mat_sd, FC, logFC, tstat, pvalall, pvaladjall),
      \(mat) {
        mat |>
          as.data.frame() |>
          tibble::rownames_to_column(feature_id_colname) |>
          return()
      }
    ) |>
      purrr::reduce(dplyr::left_join, by = feature_id_colname)
  } else {
    finalres <- purrr::map(list(FC, logFC, tstat, pvalall, pvaladjall), \(mat) {
      mat |>
        as.data.frame() |>
        tibble::rownames_to_column(feature_id_colname) |>
        return()
    }) |>
      purrr::reduce(dplyr::left_join, by = feature_id_colname)
  }

  if (return_normalized_counts == TRUE) {
    finalres <- final_res |>
      dplyr::left_join(
        v$E |>
          as.data.frame() |>
          tibble::rownames_to_column(feature_id_colname),
        by = feature_id_colname
      )
  }

  message(paste0("Total number of genes included: ", nrow(finalres)))

  ### add back Anno columns and Remove row number from Feature Column
  finalres[, feature_id_colname] <- gsub(
    "_[0-9]+$",
    "",
    finalres[, feature_id_colname]
  )
  call_me_alias <- colnames(finalres)
  colnames(finalres) <- gsub("\\(|\\)", "", call_me_alias)
  ### PH: END Create DEG Table

  ### PH: START contrast summary table input:
  #                                       -design from create Design table
  #                                       -cm from Make Contrasts
  ## Output is table showing contrasts used

  # Print out sample numbers:
  #
  sampsize <- colSums(design)
  # titleval <- "Please note Sample size:"
  # titletext <- paste(names(sampsize),
  #                    sampsize,
  #                    sep = "=",
  #                    collapse = " \n ")
  # titleall <- paste(titleval, "\n", titletext, "\n\n\n")

  contrast <- colnames(cm)
  connames <- strsplit(contrast, "-")
  connames <- lapply(connames, function(x) {
    return(gsub("\\(", "", gsub("\\)", "", x)))
  })
  contrastsize <- lapply(connames, function(x) {
    return(sampsize[unlist(x)])
  })
  footnotetext <- paste(contrast, contrastsize, sep = " : ", collapse = "\n")
  footnotetext <- paste("\n\n\nContrasts:\n", footnotetext)
  ### PH: END contrast summary table

  ### PH: START Identify DEG genes input:
  #                                   -finalres from Create DEG Table
  ## Output should be table With # of DEGs per contrast with different cutoffs
  # TODO: currently these are not used anywhere downstream
  # FCpval1 <- get_gene_lists(
  #   finalres,
  #   FC,
  #   pvalall,
  #   pvaladjall,
  #   contrasts,
  #   FClimit = 1.2,
  #   pvallimit = 0.05,
  #   pval = "pval",
  #   feature_id_colname = feature_id_colname
  # )
  # FCpval2 <- get_gene_lists(
  #   finalres,
  #   FC,
  #   pvalall,
  #   pvaladjall,
  #   contrasts,
  #   FClimit = 1.2,
  #   pvallimit = 0.01,
  #   pval = "pval",
  #   feature_id_colname = feature_id_colname
  # )
  # FCadjpval1 <- get_gene_lists(
  #   finalres,
  #   FC,
  #   pvalall,
  #   pvaladjall,
  #   contrasts,
  #   FClimit = 1.2,
  #   pvallimit = 0.05,
  #   pval = "adjpval",
  #   feature_id_colname = feature_id_colname
  # )
  # FCadjpval2 <- get_gene_lists(
  #   finalres,
  #   FC,
  #   pvalall,
  #   pvaladjall,
  #   contrasts,
  #   FClimit = 1.2,
  #   pvallimit = 0.01,
  #   pval = "adjpval",
  #   feature_id_colname = feature_id_colname
  # )
  ### PH: END Identify DEG genes

  # Mean-variance Plot.
  mv_plot <- plot_mean_variance(voom_elist = v)
  print_or_save_plot(
    mv_plot,
    filename = file.path(plots_subdir, "mean-variance.png"),
    print_plots = print_plots,
    save_plots = save_plots
  )

  df_list <- contrasts |>
    purrr::map(\(contrast) {
      finalres |>
        dplyr::select(
          tidyselect::all_of(feature_id_colname),
          tidyselect::all_of(
            purrr::map(
              contrast |>
                stringr::str_split("-") |>
                unlist() |>
                paste0(., "_"),
              tidyselect::starts_with,
              vars = colnames(.)
            ) |>
              unlist()
          ),
          tidyselect::all_of(tidyselect::starts_with(contrast))
        ) |>
        dplyr::rename_with(~ gsub(paste0(contrast, "_"), "", .x)) |>
        return()
    })

  names(df_list) <- contrasts

  moo@analyses[["diff"]] <- df_list
  return(moo)
}


get_gene_lists <- function(
  finalres,
  FC,
  pvalall,
  pvaladjall,
  contrasts,
  FClimit,
  pvallimit,
  pval,
  feature_id_colname = "Gene"
) {
  upreg_genes <- list()
  downreg_genes <- list()
  for (i in seq_len(length(contrasts))) {
    if (pval == "pval") {
      upreg_genes[[i]] <- finalres |>
        dplyr::filter(
          .data[[colnames(FC)[i]]] > FClimit &
            .data[[colnames(pvalall)[i]]] < pvallimit
        ) |>
        dplyr::pull(tidyselect::all_of(feature_id_colname)) |>
        length()
      downreg_genes[[i]] <- finalres |>
        dplyr::filter(
          .data[[colnames(FC)[i]]] < -FClimit &
            .data[[colnames(pvalall)[i]]] < pvallimit
        ) |>
        dplyr::pull(tidyselect::all_of(feature_id_colname)) |>
        length()
    } else {
      upreg_genes[[i]] <- finalres |>
        dplyr::filter(
          .data[[colnames(FC)[i]]] > FClimit &
            .data[[colnames(pvaladjall)[i]]] < pvallimit
        ) |>
        dplyr::pull(tidyselect::all_of(feature_id_colname)) |>
        length()
      downreg_genes[[i]] <- finalres |>
        dplyr::filter(
          .data[[colnames(FC)[i]]] < -FClimit &
            .data[[colnames(pvaladjall)[i]]] < pvallimit
        ) |>
        dplyr::pull(tidyselect::all_of(feature_id_colname)) |>
        length()
    }
  }
  names(upreg_genes) <- contrasts
  names(downreg_genes) <- contrasts
  allreggenes <- rbind(unlist(upreg_genes), unlist(downreg_genes))
  rownames(allreggenes) <- c(
    paste0("upreg>", FClimit, ", ", pval, "<", pvallimit),
    paste0("downreg<-", FClimit, ", ", pval, "<", pvallimit)
  )
  return(allreggenes)
}


plot_mean_variance <- function(voom_elist) {
  x <- y <- NULL
  v <- voom_elist
  sx <- v$voom.xy$x
  sy <- v$voom.xy$y
  xyplot <- as.data.frame(cbind(sx, sy))
  voomline <- as.data.frame(cbind(x = v$voom.line$x, y = v$voom.line$y))

  g <- ggplot2::ggplot() +
    ggplot2::geom_point(data = xyplot, ggplot2::aes(x = sx, y = sy), size = 1) +
    ggplot2::theme_bw() +
    ggplot2::geom_smooth(
      data = voomline,
      ggplot2::aes(x = x, y = y),
      color = "red"
    ) +
    ggplot2::ggtitle("voom: Mean-variance trend") +
    ggplot2::xlab(v$voom.xy$xlab) +
    ggplot2::ylab(v$voom.xy$ylab) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
      )
    )
  return(g)
}

#' Filter features from differential analysis based on statistical significance
#'
#' Outputs dataset of significant genes from DEG table; filters genes based on statistical significance (p-value or
#' adjusted p-value) and change (fold change, log2 fold change, or t-statistic); in addition allows for selection of DEG
#' estimates and for sub-setting of contrasts and groups included in the output gene list.
#'
#' @inheritParams option_params
#' @inheritParams filter_counts
#' @param significance_column Column name for significance, e.g. `"pval"` or `"pvaladj"` (default)
#' @param significance_cutoff Features will only be kept if their `significance_column` is less then this cutoff
#'   threshold
#' @param change_column Column name for change, e.g. `"logFC"` (default)
#' @param change_cutoff Features will only be kept if the absolute value of their `change_column` is greater than or
#'   equal to this cutoff threshold
#' @param filtering_mode Accepted values: `"any"` or `"all"` to include features that meet the criteria in _any_
#'   contrast or in _all_ contrasts
#' @param include_estimates Column names of estimates to include. Default: `c("FC", "logFC", "tstat", "pval",
#'   "adjpval")`
#' @param round_estimates Whether to round estimates. Default: `TRUE`
#' @param rounding_decimal_for_percent_cells Decimal place to use when rounding Percent cells
#' @param contrast_filter Whether to filter `contrasts` in or our of analysis. If `"keep"`, only the contrast names
#'   listed in `contrasts` will be included. If `"remove`, the contrast names listed by `contrasts` will be removed. If
#'   `"none"`, all contrasts in the dataset are used. Options: `"keep"`, `"remove"`, or `"none"`
#' @param contrasts Contrast names to filter by `contrast_filter`. If `contrast_filter` is `"none"`, this parameter has
#'   no effect.
#' @param groups Group names to filter by `groups_filter`. If `groups_filter` is `"none"`, this parameter has no effect.
#'   Options: `"keep"`, `"remove"`, or `"none"`
#' @param groups_filter Whether to filter `groups` in or out of analysis. If `"keep"`, only the group names listed in
#'   `groups` will be included. If `"remove"`, the group names listed by `groups` will be removed. If `"none"`, all
#'   groups in the dataset are used.
#' @param label_font_size Font size for labels in the plot (default: 6)
#' @param label_distance Distance of labels from the bars (default: 1)
#' @param y_axis_expansion Expansion of the y-axis (default: 0.08)
#' @param fill_colors Fill colors for the bars (default: c("steelblue1", "whitesmoke"))
#' @param pie_chart_in_3d Whether to draw pie charts in 3D (default: TRUE)
#' @param bar_width Width of the bars (default: 0.4)
#' @param draw_bar_border Whether to draw borders around bars (default: TRUE)
#' @param plot_type "bar" or "pie"
#' @param plot_titles_fontsize Font size for plot titles (default: 12)
#' @param plots_subdir subdirectory in where plots will be saved if `save_plots` is `TRUE`
#'
#' @family moo methods
#'
#' @export
#'
#' @examples
#' moo <- multiOmicDataSet(
#'   sample_metadata = as.data.frame(nidap_sample_metadata),
#'   anno_dat = data.frame(),
#'   counts_lst = list(
#'     "raw" = as.data.frame(nidap_raw_counts),
#'     "clean" = as.data.frame(nidap_clean_raw_counts),
#'     "filt" = as.data.frame(nidap_filtered_counts)
#'   )
#' ) |>
#'   diff_counts(
#'     count_type = "filt",
#'     sub_count_type = NULL,
#'     sample_id_colname = "Sample",
#'     feature_id_colname = "Gene",
#'     covariates_colnames = c("Group", "Batch"),
#'     contrast_colname = c("Group"),
#'     contrasts = c("B-A", "C-A", "B-C"),
#'     voom_normalization_method = "quantile",
#'   ) |>
#'   filter_diff()
#' head(moo@analyses$diff_filt)
filter_diff <- function(
  moo,
  feature_id_colname = NULL,
  significance_column = "adjpval",
  significance_cutoff = 0.05,
  change_column = "logFC",
  change_cutoff = 1,
  filtering_mode = "any",
  include_estimates = c("FC", "logFC", "tstat", "pval", "adjpval"),
  round_estimates = TRUE,
  rounding_decimal_for_percent_cells = 0,
  contrast_filter = "none",
  contrasts = c(),
  groups = c(),
  groups_filter = "none",
  label_font_size = 6,
  label_distance = 1,
  y_axis_expansion = 0.08,
  fill_colors = c("steelblue1", "whitesmoke"),
  pie_chart_in_3d = TRUE,
  bar_width = 0.4,
  draw_bar_border = TRUE,
  plot_type = "bar",
  plot_titles_fontsize = 12,
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = file.path("diff", "filt")
) {
  Count <- Count_format <- L1 <- Label <- Percent <- Significant <- Var1 <- Var2 <- value <- NULL

  # from NIDAP DEG_Gene_List template - filters DEG table
  diff_dat <- moo@analyses$diff |>
    join_dfs_wide() |>
    as.data.frame()

  message("* filtering differential features")

  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(diff_dat)[1]
  }
  if (!(filtering_mode %in% c("any", "all"))) {
    stop(glue::glue("filtering_mode not recognized: {filtering_mode}"))
  }
  if (!(plot_type %in% c("bar", "pie"))) {
    stop(glue::glue("plot_type not recognized: {plot_type}"))
  }
  if (!(contrast_filter %in% c("keep", "remove", "none"))) {
    stop(glue::glue("contrast_filter not recognized: {contrast_filter}"))
  }
  # ensure these are vectors, not lists. needed for reading args from JSON
  include_estimates <- include_estimates |> unlist()
  contrasts <- contrasts |> unlist()
  groups <- groups |> unlist()
  fill_colors <- fill_colors |> unlist()

  # If include_estimates param is empty, then fill it with default values.
  if (length(include_estimates) == 0) {
    include_estimates <- c("FC", "logFC", "tstat", "pval", "adjpval")
  }
  ## select DEG stat columns
  estimates <- paste0("_", include_estimates)
  signif <- paste0("_", significance_column)
  change <- paste0("_", change_column)
  diff_dat <- diff_dat |>
    dplyr::select(
      tidyselect::all_of(feature_id_colname),
      tidyselect::ends_with(c(estimates, signif, change))
    )

  contrasts_name <- diff_dat |>
    dplyr::select(tidyselect::ends_with(signif)) |>
    colnames()
  contrasts_name <- unlist(strsplit(contrasts_name, signif))
  if (contrast_filter == "keep") {
    contrasts_name <- intersect(contrasts_name, contrasts)
  } else if (contrast_filter == "remove") {
    contrasts_name <- setdiff(contrasts_name, contrasts)
  }
  contrasts_name <- paste0(contrasts_name, "_")

  groups_name <- diff_dat |>
    dplyr::select(tidyselect::ends_with(c("_mean", "_sd"))) |>
    colnames()
  groups_name <- unique(gsub("_mean|_sd", "", groups_name))
  if (groups_filter == "keep") {
    groups_name <- intersect(groups_name, groups)
  } else if (contrast_filter == "remove") {
    groups_name <- setdiff(groups_name, groups)
  }
  groups_name <- paste0(groups_name, "_")

  ### PH: END Set parameters

  ### PH: START Subset DEG table

  diff_dat <- diff_dat |>
    dplyr::select(
      tidyselect::all_of(feature_id_colname),
      tidyselect::starts_with(c(groups_name, contrasts_name))
    )

  ## select filter variables
  datsignif <- diff_dat |>
    dplyr::select(
      tidyselect::all_of(feature_id_colname),
      tidyselect::ends_with(signif)
    ) |>
    tibble::column_to_rownames(feature_id_colname)
  datchange <- diff_dat |>
    dplyr::select(
      tidyselect::all_of(feature_id_colname),
      tidyselect::ends_with(change)
    ) |>
    tibble::column_to_rownames(feature_id_colname)
  genes <- diff_dat[, feature_id_colname]

  ## filter genes
  significant <- datsignif < significance_cutoff
  changed <- abs(datchange) >= change_cutoff
  if (filtering_mode == "any") {
    selgenes <- apply(significant & changed, 1, any)
    select_genes <- genes[selgenes]
  } else if (filtering_mode == "all") {
    selgenes <- apply(significant & changed, 1, all)
    select_genes <- genes[selgenes]
  }

  # stop if 0 genes selected with the selection criteria
  if (length(select_genes) == 0) {
    stop(glue::glue(
      "ERROR: Selection criteria selected no genes - change stringency of the Significance cutoff",
      " and/or Change cutoff parameters"
    ))
  }
  message(
    glue::glue(
      "Total number of genes selected with {significance_column} < {significance_cutoff}",
      " and \u007c {change_column} \u007c \u2265 {change_cutoff} is sum(selgenes)"
    )
  )

  ## .output dataset
  out <- diff_dat |> dplyr::filter(get(feature_id_colname) %in% select_genes)
  if (round_estimates) {
    out <- out |> dplyr::mutate_if(is.numeric, ~ signif(., 3))
  }

  ### PH: END Subset DEG table

  ### PH: START Create DEG summary Barplot
  ## do plot
  significant <- apply(datsignif, 2, function(x) x <= significance_cutoff)
  changed <- apply(datchange, 2, function(x) abs(x) >= change_cutoff)
  dd <- significant & changed
  if (draw_bar_border) {
    bar_border <- "black"
  } else {
    bar_border <- NA
  }

  ## If fill_colors is blank, then
  ## give it default values.
  if (length(fill_colors) == 0) {
    fill_colors <- c("steelblue1", "whitesmoke")
  }

  if (filtering_mode == "any") {
    say_contrast <- paste(colnames(dd), collapse = " | ")
    say_contrast <- gsub("_pval|_adjpval", "", say_contrast)

    Var2df <- reshape2::melt(apply(dd, 2, table))
    if ("L1" %in% names(Var2df)) {
      Var2df <- Var2df |>
        dplyr::rename(Var2 = L1)
    }

    tab <- Var2df |>
      dplyr::mutate(Significant = ifelse(Var1, "TRUE", "FALSE")) |>
      dplyr::mutate(
        Significant = factor(Significant, levels = c("TRUE", "FALSE")),
        Count = value,
        Count_format = format(round(value, 1), nsmall = 0, big.mark = ",")
      ) |>
      dplyr::mutate(Var2 = gsub("_pval|_adjpval", "", Var2)) |>
      dplyr::group_by(Var2) |>
      dplyr::mutate(
        Percent = round(
          Count / sum(Count) * 100,
          rounding_decimal_for_percent_cells
        )
      ) |>
      dplyr::mutate(Label = sprintf("%s (%g%%)", Count_format, Percent))

    pp <- ggplot2::ggplot(
      tab,
      ggplot2::aes(
        x = "",
        y = Count,
        labels = Significant,
        fill = Significant
      )
    ) +
      ggplot2::geom_col(
        width = bar_width,
        position = "dodge",
        col = bar_border
      ) +
      ggplot2::facet_wrap(~Var2) +
      ggplot2::scale_fill_manual(values = fill_colors) +
      ggplot2::theme_bw(base_size = 20) +
      ggplot2::xlab("Contrast") +
      ggplot2::ylab("Number of Genes") +
      ggplot2::geom_text(
        ggplot2::aes(label = Label),
        color = c("black"),
        size = label_font_size,
        position = ggplot2::position_dodge(width = bar_width),
        vjust = -label_distance
      ) +
      ggplot2::theme(
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank()
      ) +
      ggplot2::ggtitle(
        sprintf(
          "%s<%g & |%s|>%g %s",
          significance_column,
          significance_cutoff,
          change_column,
          change_cutoff,
          filtering_mode
        )
      ) +
      ggplot2::theme(
        legend.key.size = ggplot2::unit(3, "line"),
        legend.position = "top",
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(size = plot_titles_fontsize)
      ) +
      ggplot2::xlab("") +
      ggplot2::scale_y_continuous(name = "", expand = c(y_axis_expansion, 0))
    print_or_save_plot(
      pp,
      filename = file.path(plots_subdir, glue::glue("{plot_type}chart.png")),
      print_plots = print_plots,
      save_plots = save_plots
    )
  } else if (filtering_mode == "all") {
    say_contrast <- paste(colnames(dd), collapse = " & ")
    say_contrast <- gsub("_pval|_adjpval", "", say_contrast)
    dd <- apply(dd, 1, function(x) all(x == TRUE))

    if (plot_type == "bar") {
      dd <- data.frame(dd)
      colnames(dd) <- say_contrast
      tab <- reshape2::melt(apply(dd, 2, table)) |>
        dplyr::mutate(Significant = ifelse(Var1, "TRUE", "FALSE")) |>
        dplyr::mutate(
          Significant = factor(Significant, levels = c("TRUE", "FALSE")),
          Count = value,
          Count_format = format(round(value, 1), nsmall = 0, big.mark = ",")
        ) |>
        dplyr::group_by(Var2) |>
        dplyr::mutate(
          Percent = round(
            Count / sum(Count) * 100,
            rounding_decimal_for_percent_cells
          )
        ) |>
        dplyr::mutate(Label = sprintf("%s (%g%%)", Count_format, Percent))

      pp <- ggplot2::ggplot(
        tab,
        ggplot2::aes(
          x = "",
          y = Count,
          labels = Significant,
          fill = Significant
        )
      ) +
        ggplot2::geom_col(
          width = bar_width,
          position = "dodge",
          col = bar_border
        ) +
        ggplot2::facet_wrap(~Var2) +
        ggplot2::scale_fill_manual(values = fill_colors) +
        ggplot2::theme_bw(base_size = 20) +
        ggplot2::xlab("Contrast") +
        ggplot2::ylab("Number of Genes") +
        ggplot2::geom_text(
          ggplot2::aes(label = Label),
          color = c("black"),
          size = label_font_size,
          position = ggplot2::position_dodge(width = bar_width),
          vjust = -label_distance
        ) +
        ggplot2::ggtitle(
          sprintf(
            "%s<%g & |%s|>%g %s",
            significance_column,
            significance_cutoff,
            change_column,
            change_cutoff,
            filtering_mode
          )
        ) +
        ggplot2::theme(
          axis.ticks.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          legend.key.size = ggplot2::unit(3, "line"),
          legend.position = "top",
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank(),
          strip.background = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(size = plot_titles_fontsize)
        ) +
        ggplot2::xlab("") +
        ggplot2::scale_y_continuous(
          name = "",
          expand = c(y_axis_expansion, 0)
        )
      print_or_save_plot(
        pp,
        filename = file.path(plots_subdir, glue::glue("{plot_type}chart.png")),
        print_plots = print_plots,
        save_plots = save_plots
      )
      ### PH: END Create DEG summary Barplot
    } else if (plot_type == "pie") {
      ### PH: START Create DEG summary PieChart
      abort_packages_not_installed("plotrix")
      N <- c(sum(dd), length(dd) - sum(dd))
      Nk <- format(round(as.numeric(N), 1), nsmall = 0, big.mark = ",")
      P <- round(N / sum(N) * 100, rounding_decimal_for_percent_cells)
      if (label_font_size > 0) {
        labs <- c(
          sprintf("Significant\n%s (%g%%)", Nk[1], P[1]),
          sprintf("Non-Significant\n%s (%g%%)", Nk[2], P[2])
        )
      } else {
        labs <- NULL
      }
      # TODO: how to print_or_save base R plot?
      if (pie_chart_in_3d) {
        plotrix::pie3D(
          N,
          radius = 0.8,
          height = 0.06,
          col = fill_colors,
          theta = 0.9,
          start = 0,
          explode = 0,
          labels = labs,
          labelcex = label_font_size,
          shade = 0.7,
          sector.order = 1:2,
          border = FALSE
        )
        graphics::title(
          main = sprintf(
            "%s<%g & |%s|>%g %s: %s",
            significance_column,
            significance_cutoff,
            change_column,
            change_cutoff,
            filtering_mode,
            say_contrast
          ),
          cex.main = plot_titles_fontsize / 3,
          line = -2
        )
      } else {
        labs <- gsub("\n", ": ", labs)
        plotrix::pie3D(
          N,
          radius = 0.8,
          height = 0.06,
          col = fill_colors,
          theta = 0.9,
          start = 45,
          explode = 0,
          labels = labs,
          labelcex = label_font_size,
          shade = 0.7,
          sector.order = 1:2,
          border = NULL
        )
        graphics::title(
          main = sprintf(
            "%s<%g & |%s|>%g %s: %s",
            significance_column,
            significance_cutoff,
            change_column,
            change_cutoff,
            filtering_mode,
            say_contrast
          ),
          cex.main = plot_titles_fontsize / 3,
          line = -2
        )
      }
    }
    ### PH: END Create DEG summary PieChart
  }
  moo@analyses$diff_filt <- out
  return(moo)
}
