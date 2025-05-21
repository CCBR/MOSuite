#' Differential expression analysis
#'
#' @inheritParams filter_counts
#' @inheritParams batch_correct_counts
#' @inheritParams normalize_counts
#' @inheritParams option_params
#'
#' @param sub_count_type if `count_type` is a list, specify the sub count type within the list. (Default: `NULL`)
#' @param covariates_colnames The column name(s) from the sample metadata
#'   containing variable(s) of interest, such as phenotype.
#'   Most commonly this will be the same column selected for your Groups Column.
#'   Some experimental designs may require that you add additional covariate columns here.
#' @param contrast_colname The column in the metadata that contains the group variables you wish to find differential expression between. Up to 2 columns (2-factor analysis) can be used.
#' @param contrasts Specify each contrast in the format group1-group2, e.g. treated-control
#' @param covariates_colnames Columns to be used as covariates in linear modeling. Must include column from "Contrast Variable". Most commonly your covariate will be group and batch (if you have different batches in your data).
#' @param return_mean_and_sd if TRUE, return Mean and Standard Deviation of groups in addition to DEG estimates for contrast(s)
#' @param return_normalized_counts if TRUE, return normalized counts for samples included in the limma model
#'
#' @returns `multiOmicDataSet` with `diff` added to the `analyses` slot (i.e. `moo@analyses$diff`)
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
#' ) %>%
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
diff_counts <- function(moo,
                        count_type = "filt",
                        sub_count_type = NULL,
                        sample_id_colname = NULL,
                        feature_id_colname = NULL,
                        samples_to_include = NULL,
                        covariates_colnames = c("Group", "Batch"), # TODO better defaults for covariates & contrasts
                        contrast_colname = c("Group"),
                        contrasts = c("B-A", "C-A", "B-C"),
                        input_in_log_counts = FALSE,
                        return_mean_and_sd = FALSE,
                        return_normalized_counts = TRUE,
                        voom_normalization_method = "quantile",
                        print_plots = options::opt("print_plots"),
                        save_plots = options::opt("save_plots"),
                        plots_subdir = "diff") {
  Sample <- group <- y <-
    sample_metadata <- moo@sample_meta
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
  message(glue::glue("* differential counts"))
  # TODO support tibbles
  counts_dat %<>% as.data.frame()


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
  counts_dat[, feature_id_colname] <- paste0(counts_dat[, feature_id_colname], "_", 1:nrow(counts_dat))

  df.m <- counts_dat
  gene_names <- NULL
  gene_names$GeneID <- counts_dat[, feature_id_colname]
  ### PH: END Check Rownames - from Filtering + Normalization Template

  ### PH: START Input Data Validation - from Filtering + Normalization Template
  ### This code block does input data validation
  # Remove samples that are not in the contrast groups:
  groups <- unique(unlist(strsplit(contrasts, "-")))
  sample_metadata <- sample_metadata %>% dplyr::filter(.data[[contrast_colname]] %in% groups)
  df.m %<>% dplyr::select(tidyr::all_of(c(
    feature_id_colname,
    sample_metadata %>% dplyr::pull(sample_id_colname)
  )))
  ### PH: END Input Data Validation - from Filtering + Normalization Template


  ####################################
  ### Computational Functions
  ################################
  ### PH: START Create Design Formula/Table
  # Put covariates in order
  covariates_colnames <- covariates_colnames[order(covariates_colnames != contrast_colname)]

  # TODO: refactor - function to sub spaces with underscores
  for (ocv in covariates_colnames) {
    sample_metadata[[ocv]] <- gsub(" ", "_", sample_metadata %>% dplyr::pull(ocv))
  }
  contrasts <- gsub(" ", "_", contrasts)
  cov <- covariates_colnames[!covariates_colnames %in% contrast_colname]

  # Combine columns if 2-factor analysis
  if (length(contrast_colname) > 1) {
    sample_metadata %>% dplyr::mutate(contmerge = paste0(.data[[contrast_colname[1]]], ".", .data[[contrast_colname[2]]])) -> sample_metadata
  } else {
    sample_metadata %>% dplyr::mutate(contmerge = .data[[contrast_colname]]) -> sample_metadata
  }

  contrast_var <- factor(sample_metadata$contmerge)

  ## create Design table
  if (length(cov) > 0) {
    dm.formula <- stats::as.formula(paste("~0 +", paste(
      "contmerge", paste(cov, sep = "+", collapse = "+"),
      sep = "+"
    )))
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
    df_unlog <- df.m %>% dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ 2^.x))
    x <- edgeR::DGEList(counts = df_unlog, genes = gene_names)
  } else {
    x <- edgeR::DGEList(counts = df.m, genes = gene_names)
  }

  # TODO add this to existing norm function & document options
  if (voom_normalization_method %in% c("TMM", "TMMwzp", "RLE", "upperquartile")) {
    x <- edgeR::calcNormFactors(x, method = voom_normalization_method)
    rownames(x) <- x$genes$GeneID
    v <- limma::voom(x, design = design, normalize = "none")
  } else {
    v <- limma::voom(x,
      design = design,
      normalize = voom_normalization_method,
      save.plot = TRUE
    )
  }
  ### PH: END Limma Normalization - Same as in Normalize Counts

  ### PH: START Linear Fit and and extract df.voom table. Could be added to Limma Normalization function above with an option to run lmFit
  rownames(v$E) <- v$genes$GeneID
  as.data.frame(v$E) %>% tibble::rownames_to_column(feature_id_colname) -> df.voom
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
    stats::p.adjust(x, "BH")
  })
  colnames(pvaladjall) <- paste(colnames(fit2$coefficients), "adjpval",
    sep =
      "_"
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
    mean.df <- as.data.frame(tve) %>%
      tibble::rownames_to_column(sample_id_colname) %>%
      dplyr::left_join(sample_metadata %>% dplyr::select(tidyr::all_of(c(sample_id_colname, contrast_colname))),
        by = sample_id_colname
      ) %>%
      dplyr::rename(group = tidyr::all_of(contrast_colname)) %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(dplyr::across(dplyr::where(is.numeric), ~ base::mean(.x))) %>%
      as.data.frame()
    mat_mean <- mean.df[, -c(1, 2)] %>%
      as.matrix() %>%
      t()
    colnames(mat_mean) <- mean.df[, 1]
    colnames(mat_mean) <- paste(colnames(mat_mean), "mean", sep = "_")
    colnames(mat_mean) <- gsub("\\.", "_", colnames(mat_mean))
    # mat_mean %<>% as.data.frame() %>% tibble::rownames_to_column(feature_id_colname)

    sd.df <- as.data.frame(tve) %>%
      tibble::rownames_to_column(sample_id_colname) %>%
      dplyr::left_join(sample_metadata %>% dplyr::select(tidyr::all_of(c(sample_id_colname, contrast_colname))),
        by = sample_id_colname
      ) %>%
      dplyr::rename(group = tidyr::all_of(contrast_colname)) %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(dplyr::across(dplyr::where(is.numeric), ~ stats::sd(.x))) %>%
      as.data.frame()
    mat_sd <- sd.df[, -c(1, 2)] %>%
      as.matrix() %>%
      t()
    colnames(mat_sd) <- sd.df[, 1]
    colnames(mat_sd) <- paste(colnames(mat_sd), "sd", sep = "_")
    colnames(mat_sd) <- gsub("\\.", "_", colnames(mat_sd))
    # mat_sd %<>% as.data.frame() %>% tibble::rownames_to_column(feature_id_colname)

    finalres <- purrr::map(list(mat_mean, mat_sd, FC, logFC, tstat, pvalall, pvaladjall), \(mat) {
      mat %>%
        as.data.frame() %>%
        tibble::rownames_to_column(feature_id_colname)
    }) %>%
      purrr::reduce(dplyr::left_join, by = feature_id_colname)
  } else {
    finalres <- purrr::map(list(FC, logFC, tstat, pvalall, pvaladjall), \(mat) {
      mat %>%
        as.data.frame() %>%
        tibble::rownames_to_column(feature_id_colname)
    }) %>%
      purrr::reduce(dplyr::left_join, by = feature_id_colname)
  }

  if (return_normalized_counts == TRUE) {
    finalres %<>% dplyr::left_join(v$E %>% as.data.frame() %>% tibble::rownames_to_column(feature_id_colname),
      by = feature_id_colname
    )
  }

  message(paste0("Total number of genes included: ", nrow(finalres)))

  ### add back Anno columns and Remove row number from Feature Column
  finalres[, feature_id_colname] <- gsub("_[0-9]+$", "", finalres[, feature_id_colname])
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
  titleval <- "Please note Sample size:"
  titletext <- paste(names(sampsize),
    sampsize,
    sep = "=",
    collapse = " \n "
  )
  titleall <- paste(titleval, "\n", titletext, "\n\n\n")

  contrast <- colnames(cm)
  connames <- strsplit(contrast, "-")
  connames <- lapply(connames, function(x) {
    gsub("\\(", "", gsub("\\)", "", x))
  })
  contrastsize <- lapply(connames, function(x) {
    sampsize[unlist(x)]
  })
  footnotetext <- paste(contrast, contrastsize, sep = " : ", collapse = "\n")
  footnotetext <- paste("\n\n\nContrasts:\n", footnotetext)
  ### PH: END contrast summary table

  ### PH: START Identify DEG genes input:
  #                                   -finalres from Create DEG Table
  ## Output should be table With # of DEGs per contrast with different cutoffs
  # TODO: currently these are not used anywhere downstream
  FCpval1 <- get_gene_lists(finalres, FC, pvalall, pvaladjall, contrasts, FClimit = 1.2, pvallimit = 0.05, pval = "pval", feature_id_colname = feature_id_colname)
  FCpval2 <- get_gene_lists(finalres, FC, pvalall, pvaladjall, contrasts, FClimit = 1.2, pvallimit = 0.01, pval = "pval", feature_id_colname = feature_id_colname)
  FCadjpval1 <- get_gene_lists(finalres, FC, pvalall, pvaladjall, contrasts, FClimit = 1.2, pvallimit = 0.05, pval = "adjpval", feature_id_colname = feature_id_colname)
  FCadjpval2 <- get_gene_lists(finalres, FC, pvalall, pvaladjall, contrasts, FClimit = 1.2, pvallimit = 0.01, pval = "adjpval", feature_id_colname = feature_id_colname)
  ### PH: END Identify DEG genes

  # Mean-variance Plot.
  mv_plot <- plot_mean_variance(voom_elist = v)
  print_or_save_plot(mv_plot,
    filename = file.path(plots_subdir, "mean-variance.png"),
    print_plots = print_plots, save_plots = save_plots
  )

  moo@analyses[["diff"]] <- finalres
  return(moo)
}


get_gene_lists <- function(finalres, FC, pvalall, pvaladjall, contrasts, FClimit, pvallimit, pval, feature_id_colname = "Gene") {
  upreg_genes <- list()
  downreg_genes <- list()
  for (i in 1:length(contrasts)) {
    if (pval == "pval") {
      finalres %>%
        dplyr::filter(.data[[colnames(FC)[i]]] > FClimit &
          .data[[colnames(pvalall)[i]]] < pvallimit) %>%
        dplyr::pull(tidyselect::all_of(feature_id_colname)) %>%
        length() -> upreg_genes[[i]]
      finalres %>%
        dplyr::filter(.data[[colnames(FC)[i]]] < -FClimit &
          .data[[colnames(pvalall)[i]]] < pvallimit) %>%
        dplyr::pull(tidyselect::all_of(feature_id_colname)) %>%
        length() -> downreg_genes[[i]]
    } else {
      finalres %>%
        dplyr::filter(.data[[colnames(FC)[i]]] > FClimit &
          .data[[colnames(pvaladjall)[i]]] < pvallimit) %>%
        dplyr::pull(tidyselect::all_of(feature_id_colname)) %>%
        length() -> upreg_genes[[i]]
      finalres %>%
        dplyr::filter(.data[[colnames(FC)[i]]] < -FClimit &
          .data[[colnames(pvaladjall)[i]]] < pvallimit) %>%
        dplyr::pull(tidyselect::all_of(feature_id_colname)) %>%
        length() -> downreg_genes[[i]]
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
    ggplot2::geom_smooth(data = voomline, ggplot2::aes(x = x, y = y), color = "red") +
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
