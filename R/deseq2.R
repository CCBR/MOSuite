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
run_deseq2 <- S7::new_generic("run_deseq2", "moo", function(moo, design, ...) {
  return(S7::S7_dispatch())
})

S7::method(run_deseq2, multiOmicDataSet) <- function(
  moo,
  design,
  feature_id_colname = "gene_id",
  min_count = 10,
  ...
) {
  if (is.null(moo@counts$filt)) {
    stop(
      "moo must contain filtered counts for DESeq2. Hint: Did you forget to run filter_counts()?"
    )
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = moo@counts$filt %>%
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), round)) %>% # DESeq2 requires integer counts
      counts_dat_to_matrix(feature_id_colname = feature_id_colname),
    colData = moo@sample_meta,
    design = design
  )
  moo@analyses$deseq2_ds <- DESeq2::DESeq(dds, ...)
  moo@analyses$deseq2_results <- DESeq2::results(moo@analyses$deseq2_ds)
  return(moo)
}
