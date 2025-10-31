#' Run DESeq2 on a multiOmicDataSet
#'
#' @param moo multiOmicDataSet object
#' @param design   model formula for experimental design. Columns must exist in `meta_dat`.
#' @param ...      remaining variables are forwarded to `DESeq2::DESeq()`.
#'
#' @return multiOmicDataSet object with DESeq2 slot filled
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' moo <- create_multiOmicDataSet_from_files(
#'   system.file("extdata", "sample_metadata.tsv.gz",
#'     package = "MOSuite"
#'   ),
#'   system.file("extdata",
#'     "RSEM.genes.expected_count.all_samples.txt.gz",
#'     package = "MOSuite"
#'   )
#' ) %>% filter_counts()
#' moo <- run_deseq2(moo, ~condition)
#' }
#' @family moo methods
run_deseq <- S7::new_generic("run_deseq2", "moo", function(moo, design, ...) {
  return(S7::S7_dispatch())
})

S7::method(run_deseq, multiOmicDataSet) <- function(moo,
                                                    design,
                                                    count_type = "filt",
                                                    sub_count_type = NULL,
                                                    feature_id_colname = NULL,
                                                    min_count = 10, ...) {
  sample_metadata <- moo@sample_meta
  message(glue::glue("* differential counts with DESeq2"))
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
  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(counts_dat)[1]
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_dat %>%
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), round)) %>% # DESeq2 requires integer counts
      counts_dat_to_matrix(feature_id_colname = feature_id_colname),
    colData = sample_metadata,
    design = design
  ) %>%
    DESeq2::DESeq(...)

  moo@analyses$deseq$dds <- dds
  moo@analyses$deseq$colData <- DESeq2::colData(dds)
  moo@analyses$deseq$res <- DESeq2::results(dds)
  moo@analyses$deseq$vsd <- DESeq2::vst(dds, blind = FALSE)
  return(moo)
}
