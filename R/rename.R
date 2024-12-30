#' Rename samples
#'
#' TODO this should happen right at the beginning of the template?
#'
#' TODO accept new names for samples in sample metadata spreadsheet
#'
#'
#' @param dat data frame
#' @param samples_to_rename_manually TODO use sample metadata spreadsheet custom column. Need to document the format of this object.
#'
#' @return data frame with samples renamed
#' @keywords internal
#'
rename_samples <- function(dat, samples_to_rename_manually) {
  replacements <- samples_to_rename_manually

  if (!is.null(replacements) && length(replacements) > 0 && nchar(replacements) > 0) {
    # TODO: refactor with dplyr::rename for simplicity
    for (x in replacements) {
      old <- strsplit(x, ": ?")[[1]][1]
      new <- strsplit(x, ": ?")[[1]][2]
      colnames(dat)[colnames(dat) %in% old] <- new
    }
  }
  return(dat)
}
