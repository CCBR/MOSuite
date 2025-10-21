#' multiOmicDataSet class
#'
#' @param sample_metadata sample metadata as a data frame or tibble. The first column is assumed to contain the sample
#'   IDs which must correspond to column names in the raw counts.
#' @param anno_dat data frame of feature annotations, such as gene symbols or any other information about the features
#'   in `counts_lst`.
#' @param counts_lst named list of data frames containing counts, e.g. expected feature counts from RSEM. Each data
#'   frame is expected to contain a `feature_id` column as the first column, and all remaining columns are sample IDs in
#'   the `sample_meta`.
#' @param analyses_lst named list of analysis results, e.g. DESeq results object
#' @export
#'
#' @family moo constructors
multiOmicDataSet <- S7::new_class(
  "multiOmicDataSet",
  properties = list(
    sample_meta = S7::class_data.frame,
    annotation = S7::class_data.frame,
    counts = S7::class_list,
    # list of data frames
    analyses = S7::class_list
  ),
  constructor = function(sample_metadata,
                         anno_dat,
                         counts_lst,
                         analyses_lst = list()) {
    if (!("colors" %in% names(analyses_lst))) {
      analyses_lst[["colors"]] <- get_colors_lst(sample_metadata)
    }
    return(S7::new_object(
      S7::S7_object(),
      sample_meta = sample_metadata,
      annotation = anno_dat,
      counts = counts_lst,
      analyses = analyses_lst
    ))
  },
  validator = function(self) {
    # counts must only contain approved names
    approved_counts <- c("raw", "clean", "cpm", "filt", "norm", "batch")
    if (!all(names(self@counts) %in% approved_counts)) {
      stop(
        glue::glue(
          "counts can only contain these names:\n\t{paste(approved_counts, collapse = ', ')}"
        )
      )
    }

    if (!("raw" %in% names(self@counts))) {
      stop("counts must contain at least 'raw' counts")
    }
    meta_sample_colnames <- self@sample_meta %>% dplyr::pull(1)
    feature_sample_colnames <- self@counts$raw %>%
      dplyr::select(-1) %>%
      colnames()

    # all sample IDs in sample_meta must also be in raw counts, & vice versa
    error_msg <- ""
    not_in_meta <- setdiff(meta_sample_colnames, feature_sample_colnames)
    if (length(not_in_meta) > 0) {
      error_msg <- glue::glue(
        "Not all columns after the first column in the raw counts data are sample IDs in the sample metadata:\n\t",
        "{glue::glue_collapse(not_in_meta, sep = ', ')}"
      )
    }
    not_in_counts <- setdiff(feature_sample_colnames, meta_sample_colnames)
    if (length(not_in_counts) > 0) {
      error_msg <- glue::glue(
        "Not all sample IDs in the sample metadata are in the raw count data:\n\t",
        "{glue::glue_collapse(not_in_counts, sep = ', ')}"
      )
    }
    if (nchar(error_msg) > 0) {
      stop(error_msg)
    }

    # sample IDs must be in the same order
    if (!all(feature_sample_colnames == meta_sample_colnames)) {
      stop(glue::glue(
        "The sample IDs in the sample metadata do not equal the columns in the raw count data. ",
        "Sample IDs must be in the same order."
      ))
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
#' @param sample_id_colname name of the column in `sample_metadata` that contains the sample IDs. (Default: `NULL` -
#'   first column in the sample metadata will be used.)
#' @param feature_id_colname name of the column in `counts_dat` that contains feature/gene IDs. (Default: `NULL` - first
#'   column in the count data will be used.)
#'
#' @return [multiOmicDataSet] object
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
#'
#' sample_meta_nidap <- readr::read_csv(system.file("extdata", "nidap",
#'   "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz",
#'   package = "MOSuite"
#' ))
#' raw_counts_nidap <- readr::read_csv(system.file("extdata", "nidap", "Raw_Counts.csv.gz", package = "MOSuite"))
#' moo_nidap <- create_multiOmicDataSet_from_dataframes(sample_meta_nidap, raw_counts_nidap)
#'
#' @family moo constructors
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
    stop(
      glue::glue(
        "Not all sample IDs in the sample metadata are in the count data. Samples missing in count data:\n\t",
        glue::glue_collapse(meta_sample_colnames[!(meta_sample_colnames %in% colnames(counts_dat))], sep = ", ")
      )
    )
  }

  # create anno_dat out of excess columns in count dat
  anno_dat <- counts_dat %>%
    dplyr::select(-tidyselect::all_of(meta_sample_colnames))
  counts_dat <- counts_dat %>%
    dplyr::select(
      !!rlang::sym(feature_id_colname),
      tidyselect::all_of(meta_sample_colnames)
    )

  counts <- list()
  counts[[count_type]] <- counts_dat

  return(multiOmicDataSet(sample_metadata, anno_dat, counts))
}

#' Construct a multiOmicDataSet object from tsv files.
#'
#' @inheritParams multiOmicDataSet
#' @inheritParams create_multiOmicDataSet_from_dataframes
#' @param sample_meta_filepath path to text file with sample IDs and metadata for differential analysis.
#' @param feature_counts_filepath path to text file of expected feature counts (e.g. gene counts from RSEM).
#' @param ... additional arguments forwarded to `readr::read_delim()`
#'
#' @return [multiOmicDataSet] object
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
#'
#' moo_nidap <- create_multiOmicDataSet_from_files(
#'   system.file("extdata", "nidap",
#'     "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz",
#'     package = "MOSuite"
#'   ),
#'   system.file("extdata", "nidap", "Raw_Counts.csv.gz", package = "MOSuite"),
#'   delim = ","
#' )
#'
#' @family moo constructors
create_multiOmicDataSet_from_files <- function(sample_meta_filepath,
                                               feature_counts_filepath,
                                               count_type = "raw",
                                               sample_id_colname = NULL,
                                               feature_id_colname = NULL,
                                               ...) {
  counts_dat <- readr::read_delim(feature_counts_filepath, ...)
  sample_metadata <- readr::read_delim(sample_meta_filepath, ...)
  return(
    create_multiOmicDataSet_from_dataframes(
      sample_metadata = sample_metadata,
      counts_dat = counts_dat,
      count_type = "raw",
      sample_id_colname = sample_id_colname,
      feature_id_colname = feature_id_colname
    )
  )
}

#' Extract count data
#'
#' @usage
#' extract_counts(moo, count_type, sub_count_type = NULL)
#'
#' @param moo multiOmicDataSet containing `count_type` & `sub_count_type` in the counts slot
#' @param count_type the type of counts to use -- must be a name in the counts slot (`moo@counts[[count_type]]`)
#' @param sub_count_type if `count_type` is a list, specify the sub count type within the list
#'   (`moo@counts[[count_type]][[sub_count_type]]`). (Default: `NULL`)
#'
#' @export
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
#' )
#'
#' moo %>%
#'   extract_counts("filt") %>%
#'   head()
#'
#' moo %>%
#'   extract_counts("norm", "voom") %>%
#'   head()
#'
extract_counts <- S7::new_generic("extract_counts", "moo", function(moo, count_type, sub_count_type = NULL) {
  return(S7::S7_dispatch())
})

S7::method(extract_counts, multiOmicDataSet) <- function(moo, count_type, sub_count_type = NULL) {
  # select correct counts matrix
  if (!(count_type %in% names(moo@counts))) {
    stop(
      glue::glue(
        "count_type {count_type} not in moo@counts. Count types: {glue::glue_collapse(names(moo@counts), sep = ', ')}"
      )
    )
  }
  counts_dat <- moo@counts[[count_type]]
  if (!is.null(sub_count_type)) {
    if (!(inherits(counts_dat, "list"))) {
      stop(
        glue::glue(
          "{count_type} counts does not contain subtypes. To use {count_type} counts, set sub_count_type to NULL"
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
  } else if (inherits(counts_dat, "list")) {
    stop(
      glue::glue(
        "{count_type} counts contains subtypes. You must set sub_count_type to extract a subtype"
      )
    )
  }
  return(counts_dat)
}

#' @name moo_counts
#'
#' @description
#'
#' The first argument can be a `multiOmicDataset` object (`moo`) or a `data.frame` containing counts.
#' For a `moo`, choose which counts slot to use with `count_type` & (optionally) `sub_count_type`.
#' For a `data.frame`, you must also set `sample_metadata`.
#' All other arguments are optional.
NULL
