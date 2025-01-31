## training data set from the NIDAP Bulk RNA-seq workflow

nidap_sample_metadata <- readr::read_csv(system.file(
  "extdata", "nidap",
  "Sample_Metadata_Bulk_RNA-seq_Training_Dataset_CCBR.csv.gz",
  package = "MOSuite"
))
usethis::use_data(nidap_sample_metadata, overwrite = TRUE)

nidap_raw_counts <- readr::read_csv(system.file(
  "extdata", "nidap", "Raw_Counts.csv.gz",
  package = "MOSuite"
))
usethis::use_data(nidap_raw_counts, overwrite = TRUE)

nidap_clean_raw_counts <- readr::read_csv(system.file(
  "extdata", "nidap", "Clean_Raw_Counts.csv.gz",
  package = "MOSuite"
))
usethis::use_data(nidap_clean_raw_counts, overwrite = TRUE)

nidap_filtered_counts <- readr::read_csv(system.file(
  "extdata", "nidap",
  "Filtered_Counts.csv.gz",
  package = "MOSuite"
))
usethis::use_data(nidap_filtered_counts, overwrite = TRUE)

nidap_norm_counts <- readr::read_csv(system.file(
  "extdata", "nidap", "Normalized_Counts.csv.gz",
  package = "MOSuite"
))
usethis::use_data(nidap_norm_counts, overwrite = TRUE)

nidap_batch_corrected_counts <- readr::read_csv(system.file("extdata", "nidap", "Batch_Corrected_Counts.csv.gz", package = "MOSuite"))
usethis::use_data(nidap_batch_corrected_counts, overwrite = TRUE)

moo <- multiOmicDataSet(
  sample_meta_dat = as.data.frame(nidap_sample_metadata),
  anno_dat = data.frame(),
  counts_lst = list(
    "raw" = as.data.frame(nidap_raw_counts),
    "norm" = list(
      "voom" = as.data.frame(nidap_norm_counts)
    )
  )
) %>%
  batch_correct_counts(
    count_type = "norm",
    sub_count_type = "voom",
    covariates_colname = "Group",
    batch_colname = "Batch",
    label_colname = "Label"
  )
nidap_batch_corrected_counts_2 <- moo@counts[["batch"]]
usethis::use_data(nidap_batch_corrected_counts_2, overwrite = TRUE)
