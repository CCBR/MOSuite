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
