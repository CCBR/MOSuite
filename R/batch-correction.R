#' Perform batch correction
#'
#' @inheritParams filter_counts
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
#'   sample_meta_dat = as.data.frame(nidap_sample_metadata),
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
#'     covariates_colname = "Group",
#'     batch_colname = "Batch",
#'     label_column = "Label"
#'   )
#'
#' head(moo@counts[["batch"]])
#'
batch_correct_counts <- function(moo,
                                 count_type = "norm", sub_count_type = "voom",
                                 sample_id_colname = NULL,
                                 feature_id_colname = NULL,
                                 samples_to_include = NULL,
                                 covariates_colnames = "Group",
                                 batch_colname = "Batch",
                                 label_column = "Label") {
  # select correct counts matrix
  if (!(count_type %in% names(moo@counts))) {
    stop(glue::glue("count_type {count_type} not in moo@counts"))
  }
  counts_dat <- moo@counts[[count_type]]
  if (!is.null(sub_count_type)) {
    if (!(inherits(counts_dat, "list"))) {
      stop(glue::glue("{count_type} counts is not a named list. To use {count_type} counts, set sub_count_type to NULL"))
    } else if (!(sub_count_type %in% names(counts_dat))) {
      stop(glue::glue("sub_count_type {sub_count_type} is not in moo@counts[[{count_type}]]"))
    }
    counts_dat <- moo@counts[[count_type]][[sub_count_type]]
  }
  sample_metadata <- moo@sample_meta
  batch_vctr <- sample_metadata %>% dplyr::pull(batch_colname)


  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(counts_dat)[1]
  }
  if (is.null(samples_to_include)) {
    samples_to_include <- sample_metadata %>% dplyr::pull(sample_id_colname)
  }

  if (batch_colname %in% covariates_colnames) {
    stop(glue::glue("Batch column ({batch_colname}) cannot be included in covariates."))
  }
  counts_matr <- counts_dat %>%
    counts_dat_to_matrix(feature_id_colname = feature_id_colname)

  if (length(unique(batch_vctr)) <= 1) {
    combat_edata <- counts_matr
    warning(glue::glue("Batch column {batch_column} contains only 1 unique value; skipping batch correction"))
  } else {
    for (cov in covariates_colnames) {
      # TODO use mutate across
      sample_metadata[[cov]] <- as.factor(sample_metadata[[cov]])
    }
    dm.formula <- as.formula(paste("~", paste(covariates_colnames, sep = "+", collapse = "+")))
    modcombat <- model.matrix(dm.formula, data = sample_metadata)
    combat_edata <- sva::ComBat(
      counts_matr,
      batch = batch_vctr,
      mod = modcombat,
      par.prior = TRUE,
      prior.plots = FALSE
    )
  }

  combat_edata <- as.data.frame(combat_edata) %>%
    tibble::rownames_to_column(feature_id_colname)

  message(glue::glue("The total number of features in output: {nrow(combat_edata)}"))
  message(glue::glue("Samples:\n\t{colnames(combat_edata[, !colnames(combat_edata) %in% feature_id_colname])}"))
  message(glue::glue("Number of samples after batch correction: {length(samples_to_include)}"))

  moo@counts[["batch"]] <- combat_edata
  return(moo)
}
