#' multiOmicDataSet class
#'
#' @param sample_metadata sample metadata as a data frame or tibble.
#'   The first column is assumed to contain the sample IDs which must correspond to column names in the raw counts.
#' @param anno_dat data frame of feature annotations, such as gene symbols or any other information about the features in `counts_lst`.
#' @param counts_lst named list of data frames containing counts, e.g. expected feature counts from RSEM.
#'   Each data frame is expected to contain a `feature_id` column as the first column, and all remaining columns are sample IDs in the `sample_meta`.
#' @param analyses_lst named list of analysis results, e.g. DESeq results object
#' @export
multiOmicDataSet <- S7::new_class("multiOmicDataSet",
  properties = list(
    sample_meta = S7::class_data.frame,
    annotation = S7::class_data.frame,
    counts = S7::class_list, # list of data frames
    analyses = S7::class_list
  ),
  constructor = function(sample_metadata, anno_dat, counts_lst, analyses_lst = list()) {
    if (!("colors" %in% names(analyses_lst))) {
      analyses_lst[["colors"]] <- set_colors(sample_metadata)
    }
    S7::new_object(S7::S7_object(),
      sample_meta = sample_metadata,
      annotation = anno_dat,
      counts = counts_lst,
      analyses = analyses_lst
    )
  },
  validator = function(self) {
    # counts must only contain approved names
    approved_counts <- c("raw", "clean", "cpm", "filt", "norm", "batch")
    if (!all(names(self@counts) %in% approved_counts)) {
      stop(glue::glue("counts can only contain these names:\n\t{paste(approved_counts, collapse = ', ')}"))
    }

    # all sample IDs in sample_meta must also be in raw counts, & vice versa
    meta_sample_colnames <- self@sample_meta %>% dplyr::pull(1)
    feature_sample_colnames <- self@counts$raw %>%
      dplyr::select(-1) %>%
      colnames()

    error_msg <- ""
    not_in_meta <- setdiff(meta_sample_colnames, feature_sample_colnames)
    if (length(not_in_meta) > 0) {
      error_msg <- glue::glue("Not all columns after the first column in the raw counts data are sample IDs in the sample metadata:\n\t{glue::glue_collapse(not_in_meta, sep = ', ')}")
    }
    not_in_counts <- setdiff(feature_sample_colnames, meta_sample_colnames)
    if (length(not_in_counts) > 0) {
      error_msg <- glue::glue("Not all sample IDs in the sample metadata are in the raw count data:\n\t{glue::glue_collapse(not_in_counts, sep = ', ')}")
    }
    if (nchar(error_msg) > 0) {
      stop(error_msg)
    }

    # sample IDs must be in the same order
    if (!all(feature_sample_colnames == meta_sample_colnames)) {
      stop("The sample IDs in the sample metadata do not equal the columns in the raw count data. Sample IDs must be in the same order.")
    }

    # TODO any sample ID in filt or norm_cpm counts must also be in sample_meta
    # TODO counts can only contain 1 feature name column, and all other columns are sample counts
  }
)

#' Construct a multiOmicDataSet object from data frames
#'
#' @inheritParams multiOmicDataSet
#' @param counts_dat data frame of feature counts (e.g. expected feature counts from RSEM).
#' @param count_type type to assign the values of `counts_dat` to in the `counts` slot
#' @param sample_id_colname name of the column in `sample_metadata` that contains the sample IDs. (Default: `NULL` - first column in the sample metadata will be used.)
#' @param feature_id_colname name of the column in `counts_dat` that contains feature/gene IDs. (Default: `NULL` - first column in the count data will be used.)
#'
#' @return multiOmicDataSet object
#' @export
#'
#' @examples
#' sample_meta <- data.frame(
#'   sample_id = c("KO_S3", "KO_S4", "WT_S1", "WT_S2"),
#'   condition = factor(
#'     c("knockout", "knockout", "wildtype", "wildtype"),
#'     levels = c("wildtype", "knockout")
#'   )
#' )
#' moo <- create_multiOmicDataSet_from_dataframes(sample_meta, gene_counts)
#' head(moo@sample_meta)
#' head(moo@counts$raw)
#' head(moo@annotation)
create_multiOmicDataSet_from_dataframes <- function(sample_metadata,
                                                    counts_dat,
                                                    sample_id_colname = NULL,
                                                    feature_id_colname = NULL,
                                                    count_type = "raw") {
  # move sample & feature ID columns to first
  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  } else {
    sample_metadata <- sample_metadata %>%
      dplyr::relocate(!!rlang::sym(sample_id_colname))
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(counts_dat)[1]
  } else {
    counts_dat <- counts_dat %>%
      dplyr::relocate(!!rlang::sym(feature_id_colname))
  }

  meta_sample_colnames <- sample_metadata %>% dplyr::pull(sample_id_colname)
  if (!all(meta_sample_colnames %in% colnames(counts_dat))) {
    stop("Not all sample IDs in the sample metadata are in the count data")
  }
  feature_sample_colnames <- counts_dat %>%
    dplyr::select(tidyselect::all_of(meta_sample_colnames)) %>%
    colnames()

  # create anno_dat out of excess columns in count dat
  anno_dat <- counts_dat %>%
    dplyr::select(-tidyselect::all_of(meta_sample_colnames))
  counts_dat <- counts_dat %>%
    dplyr::select(!!rlang::sym(feature_id_colname), tidyselect::all_of(meta_sample_colnames))

  counts <- list()
  counts[[count_type]] <- counts_dat

  return(multiOmicDataSet(sample_metadata, anno_dat, counts))
}

#' Construct a multiOmicDataSet object from tsv files.
#'
#' @inheritParams multiOmicDataSet
#' @inheritParams create_multiOmicDataSet_from_dataframes
#' @param sample_meta_filepath path to tsv file with sample IDs and metadata for differential analysis.
#' @param feature_counts_filepath path to tsv file of expected feature counts (e.g. gene counts from RSEM).
#' @return multiOmicDataSet object
#' @export
#'
#' @examples
#' moo <- create_multiOmicDataSet_from_files(
#'   sample_meta_filepath = system.file("extdata",
#'     "sample_metadata.tsv.gz",
#'     package = "MOSuite"
#'   ),
#'   feature_counts_filepath = system.file("extdata",
#'     "RSEM.genes.expected_count.all_samples.txt.gz",
#'     package = "MOSuite"
#'   )
#' )
#' moo@counts$raw %>% head()
#' moo@sample_meta
create_multiOmicDataSet_from_files <- function(sample_meta_filepath, feature_counts_filepath,
                                               count_type = "raw",
                                               sample_id_colname = NULL,
                                               feature_id_colname = NULL) {
  counts_dat <- readr::read_tsv(feature_counts_filepath)
  sample_metadata <- readr::read_tsv(sample_meta_filepath)
  return(create_multiOmicDataSet_from_dataframes(
    sample_metadata = sample_metadata,
    counts_dat = counts_dat,
    count_type = "raw",
    sample_id_colname = sample_id_colname,
    feature_id_colname = feature_id_colname
  ))
}
