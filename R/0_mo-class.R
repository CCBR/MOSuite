# Re-exports and wrappers for MOObject::multiOmicDataSet and related functions

#' @importFrom MOObject multiOmicDataSet
#' @export
MOObject::multiOmicDataSet

#' Construct a multiOmicDataSet object from data frames
#'
#' @inheritParams MOObject::create_multiOmicDataSet_from_dataframes
#' @seealso [MOObject::create_multiOmicDataSet_from_dataframes()]
#' @returns [multiOmicDataSet] object
#' @export
#' @family moo constructors
create_multiOmicDataSet_from_dataframes <- function(
  sample_metadata,
  counts_dat,
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  count_type = "raw"
) {
  moo <- MOObject::create_multiOmicDataSet_from_dataframes(
    sample_metadata = sample_metadata,
    counts_dat = counts_dat,
    sample_id_colname = sample_id_colname,
    feature_id_colname = feature_id_colname,
    count_type = count_type
  )

  if (!("colors" %in% names(moo@analyses))) {
    moo@analyses[["colors"]] <- get_colors_lst(sample_metadata)
  }

  return(moo)
}

#' Construct a multiOmicDataSet object from text files (e.g. TSV, CSV).
#'
#' @inheritParams MOObject::create_multiOmicDataSet_from_files
#' @seealso [MOObject::create_multiOmicDataSet_from_files()]
#' @returns [multiOmicDataSet] object
#' @export
#' @family moo constructors
create_multiOmicDataSet_from_files <- function(
  sample_meta_filepath,
  feature_counts_filepath,
  count_type = "raw",
  sample_id_colname = NULL,
  feature_id_colname = NULL,
  delim = NULL,
  ...
) {
  moo <- MOObject::create_multiOmicDataSet_from_files(
    sample_meta_filepath = sample_meta_filepath,
    feature_counts_filepath = feature_counts_filepath,
    count_type = count_type,
    sample_id_colname = sample_id_colname,
    feature_id_colname = feature_id_colname,
    delim = delim,
    ...
  )

  if (!("colors" %in% names(moo@analyses))) {
    moo@analyses[["colors"]] <- get_colors_lst(moo@sample_meta)
  }

  return(moo)
}

#' Extract count data
#'
#' @inheritParams MOObject::extract_counts
#' @seealso [MOObject::extract_counts()]
#' @export
extract_counts <- MOObject::extract_counts

#' Write a multiOmicDataSet to disk as an RDS file
#'
#' @inheritParams MOObject::write_multiOmicDataSet
#' @seealso [MOObject::write_multiOmicDataSet()]
#' @return Invisibly returns `filepath`
#' @export
write_multiOmicDataSet <- MOObject::write_multiOmicDataSet

#' Read a multiOmicDataSet from disk
#'
#' @inheritParams MOObject::read_multiOmicDataSet
#' @seealso [MOObject::read_multiOmicDataSet()]
#' @return [multiOmicDataSet]
#' @export
read_multiOmicDataSet <- MOObject::read_multiOmicDataSet

#' Write multiOmicDataSet properties to disk as CSV files
#'
#' @inheritParams MOObject::write_multiOmicDataSet_properties
#' @seealso [MOObject::write_multiOmicDataSet_properties()]
#' @return Invisibly returns the `output_dir` where the files were saved
#' @export
write_multiOmicDataSet_properties <- MOObject::write_multiOmicDataSet_properties
