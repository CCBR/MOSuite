#' Validate raw counts sample metadata for filtering
#'
#' @inheritParams filter_counts
#'
#' @return sample metadata with empty cells removed and special characters replaced
#' @keywords internal
#'
#'
validate_sample_metadata <- function(counts_matrix,
                                     sample_metadata,
                                     sample_names_column = "Sample",
                                     group_column = "Group") {
  sample_metadata <- sample_metadata[match(colnames(counts_matrix), sample_metadata[[sample_names_column]]), ] # First match sample metadata to counts matrix
  sample_metadata <- sample_metadata[rowSums(is.na(sample_metadata)) != ncol(sample_metadata), ] # Remove empty rows
  sample_metadata <- sample_metadata[, colSums(is.na(sample_metadata)) == 0] # Remove empty columns
  rownames(sample_metadata) <- sample_metadata[[sample_names_column]]


  ### Remove special characters from Metadata Column. Replace with _
  sample_metadata[, group_column] <- gsub("-| |!|\\*|\\.", "_", sample_metadata[, group_column])

  return(sample_metadata)
}
