#' multiOmicDataSet class
#'
#' @param sample_meta_dat sample metadata as a data frame or tibble.
#'   The first column is assumed to contain the sample IDs which must correspond to column names in the raw counts.
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
  constructor = function(sample_meta_dat, anno_dat, counts_lst, analyses_lst = list()) {
    S7::new_object(S7::S7_object(),
      sample_meta = sample_meta_dat,
      annotation = anno_dat,
      counts = counts_lst,
      analyses = analyses_lst
    )
  },
  validator = function(self) {
    # counts must only contain approved names
    approved_counts <- c("raw", "clean", "cpm", "filt", "norm")
    if (!all(names(self@counts) %in% approved_counts)) {
      stop(glue::glue("counts can only contain these names:\n\t{paste(approved_counts, collapse = ', ')}"))
    }

    # TODO all sample IDs must be in both sample_meta and raw counts
    # TODO any sample ID in filt or norm_cpm counts must also be in sample_meta
    # TODO counts can only contain 1 feature name column, and all other columns are sample counts
  }
)

#' Construct a multiOmicDataSet object from data frames
#'
#' @inheritParams multiOmicDataSet
#' @param count_dat data frame of feature counts (e.g. expected feature counts from RSEM).
#' @param count_type type to assign the values of `count_dat` to in the `counts` slot
#' @param sample_id_colname name of the column in `sample_meta_dat` that contains the sample IDs
#' @param feature_id_colname name of the column in `count_dat` that contains feature/gene IDs
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
create_multiOmicDataSet_from_dataframes <- function(sample_meta_dat,
                                                    count_dat,
                                                    sample_id_colname = "sample_id",
                                                    feature_id_colname = "gene_id",
                                                    count_type = "raw") {
  # move sample & feature ID columns to first
  count_dat <- count_dat %>%
    dplyr::relocate(!!rlang::sym(feature_id_colname))
  sample_meta_dat <- sample_meta_dat %>%
    dplyr::relocate(!!rlang::sym(sample_id_colname))

  # make sure sample IDs are in the counts data
  meta_sample_colnames <- sample_meta_dat %>% dplyr::pull(sample_id_colname)
  if (!all(meta_sample_colnames %in% colnames(count_dat))) {
    stop("Not all sample IDs in the sample metadata are in the count data.")
  }
  feature_sample_colnames <- count_dat %>%
    dplyr::select(tidyselect::all_of(meta_sample_colnames)) %>%
    colnames()
  # sample IDs must be in the same order
  if (!all(feature_sample_colnames == meta_sample_colnames)) {
    stop("Not all columns in the count data equal the rows in the first column of the sample metadata. Sample IDs must be in the same order.")
  }

  # create anno_dat out of excess columns in count dat
  anno_dat <- count_dat %>%
    dplyr::select(-tidyselect::all_of(meta_sample_colnames))
  count_dat <- count_dat %>%
    dplyr::select(!!rlang::sym(feature_id_colname), tidyselect::all_of(meta_sample_colnames))

  counts <- list()
  counts[[count_type]] <- count_dat

  return(multiOmicDataSet(sample_meta_dat, anno_dat, counts))
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
#'     "RSEM.features.expected_count.all_samples.txt.gz",
#'     package = "MOSuite"
#'   )
#' )
#' moo@counts$raw %>% head()
#' moo@sample_meta
create_multiOmicDataSet_from_files <- function(sample_meta_filepath, feature_counts_filepath,
                                               count_type = "raw",
                                               sample_id_colname = "sample_id",
                                               feature_id_colname = "gene_id") {
  count_dat <- readr::read_tsv(feature_counts_filepath)
  sample_meta_dat <- readr::read_tsv(sample_meta_filepath)
  return(create_multiOmicDataSet_from_dataframes(
    sample_meta_dat = sample_meta_dat,
    count_dat = count_dat,
    count_type = "raw",
    sample_id_colname = sample_id_colname,
    feature_id_colname = feature_id_colname
  ))
}
