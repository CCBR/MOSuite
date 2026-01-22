#' Perform batch correction
#'
#' Perform batch correction using sva::ComBat()
#'
#' @inheritParams filter_counts
#' @inheritParams option_params
#'
#' @param sub_count_type if `count_type` is a list, specify the sub count type within the list. (Default: `"voom"`)
#' @param covariates_colnames The column name(s) from the sample metadata
#'   containing variable(s) of interest, such as phenotype.
#'   Most commonly this will be the same column selected for your Groups Column.
#'   Some experimental designs may require that you add additional covariate columns here.
#'   Do not include the `batch_colname` here.
#' @param batch_colname The column from the sample metadata containing the batch information.
#'   Samples extracted, prepared, or sequenced at separate times or using separate materials/staff/equipment
#'   may belong to different batches.
#'   Not all data sets have batches, in which case you do not need batch correction.
#'   If your data set has no batches, you can provide a batch column with the same
#'   value in every row to skip batch correction (alternatively, simply do not run this function).
#'
#' @return `multiOmicDataSet` with batch-corrected counts
#' @export
#'
#' @examples
#' moo <- multiOmicDataSet(
#'   sample_metadata = as.data.frame(nidap_sample_metadata),
#'   anno_dat = data.frame(),
#'   counts_lst = list(
#'     "raw" = as.data.frame(nidap_raw_counts),
#'     "clean" = as.data.frame(nidap_clean_raw_counts),
#'     "filt" = as.data.frame(nidap_filtered_counts),
#'     "norm" = list(
#'       "voom" = as.data.frame(nidap_norm_counts)
#'     )
#'   )
#' ) %>%
#'   batch_correct_counts(
#'     count_type = "norm",
#'     sub_count_type = "voom",
#'     covariates_colnames = "Group",
#'     batch_colname = "Batch",
#'     label_colname = "Label"
#'   )
#'
#' head(moo@counts[["batch"]])
#'
#' @family moo methods
batch_correct_counts <- function(
  moo,
  count_type = "norm",
  sub_count_type = "voom",
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  samples_to_include = NULL,
  covariates_colnames = "Group",
  batch_colname = "Batch",
  label_colname = NULL,
  colors_for_plots = NULL,
  print_plots = options::opt("print_plots"),
  save_plots = options::opt("save_plots"),
  plots_subdir = "batch"
) {
  abort_packages_not_installed("sva")
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
  # sva::ComBat() can't handle tibbles
  counts_dat <- counts_dat %>% as.data.frame()
  sample_metadata <- moo@sample_meta %>% as.data.frame()
  batch_vctr <- sample_metadata %>% dplyr::pull(batch_colname)
  message(
    glue::glue(
      "* batch-correcting {glue::glue_collapse(c(count_type, sub_count_type),sep='-')} counts"
    )
  )

  covariates_colnames <- covariates_colnames %>% unlist()

  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(counts_dat)[1]
  }
  if (is.null(samples_to_include)) {
    samples_to_include <- sample_metadata %>% dplyr::pull(sample_id_colname)
  }
  if (is.null(label_colname)) {
    label_colname <- sample_id_colname
  }

  if (batch_colname %in% covariates_colnames) {
    stop(glue::glue(
      "Batch column ({batch_colname}) cannot be included in covariates."
    ))
  }
  if (length(unique(batch_vctr)) <= 1) {
    combat_edata <- counts_dat
    warning(
      glue::glue(
        "Batch column {batch_column} contains only 1 unique value; skipping batch correction"
      )
    )
  } else {
    counts_matr <- counts_dat %>%
      counts_dat_to_matrix(feature_id_colname = feature_id_colname)
    # coerce covariate columns to factors
    sample_metadata <- sample_metadata %>%
      dplyr::mutate(dplyr::across(
        tidyselect::all_of(covariates_colnames),
        ~ as.factor(.x)
      ))
    # run batch correction
    combat_edata <- sva::ComBat(
      counts_matr,
      batch = batch_vctr,
      mod = stats::model.matrix(
        stats::as.formula(paste(
          "~",
          paste(covariates_colnames, sep = "+", collapse = "+")
        )),
        data = sample_metadata
      ),
      par.prior = TRUE,
      prior.plots = FALSE
    ) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(feature_id_colname)
  }

  if (isTRUE(print_plots) || isTRUE(save_plots)) {
    if (is.null(colors_for_plots)) {
      colors_for_plots <- moo@analyses[["colors"]][[batch_colname]]
    }
    pca_plot <- plot_pca(
      combat_edata,
      sample_metadata = sample_metadata,
      sample_id_colname = sample_id_colname,
      feature_id_colname = feature_id_colname,
      group_colname = batch_colname,
      label_colname = label_colname,
      color_values = colors_for_plots
    ) +
      ggplot2::labs(caption = "batch-corrected counts")

    hist_plot <- plot_histogram(
      combat_edata,
      sample_metadata,
      sample_id_colname = sample_id_colname,
      feature_id_colname = feature_id_colname,
      group_colname = batch_colname,
      label_colname = label_colname,
      color_values = colors_for_plots,
      color_by_group = TRUE
    ) +
      ggplot2::labs(caption = "batch-corrected counts")
    corHM_plot <- plot_corr_heatmap(
      combat_edata,
      sample_metadata = sample_metadata,
      sample_id_colname = sample_id_colname,
      feature_id_colname = feature_id_colname,
      group_colname = batch_colname,
      label_colname = label_colname,
      color_values = colors_for_plots
    ) +
      ggplot2::labs(caption = "batch-corrected counts")

    print_or_save_plot(
      pca_plot,
      filename = file.path(plots_subdir, "pca.png"),
      print_plots = print_plots,
      save_plots = save_plots
    )
    print_or_save_plot(
      hist_plot,
      filename = file.path(plots_subdir, "histogram.png"),
      print_plots = print_plots,
      save_plots = save_plots
    )
    print_or_save_plot(
      corHM_plot,
      filename = file.path(plots_subdir, "corr_heatmap.png"),
      print_plots = print_plots,
      save_plots = save_plots
    )
  }

  message(glue::glue(
    "The total number of features in output: {nrow(combat_edata)}"
  ))
  message(glue::glue(
    "Number of samples after batch correction: {ncol(combat_edata)}"
  ))

  moo@counts[["batch"]] <- combat_edata
  return(moo)
}
