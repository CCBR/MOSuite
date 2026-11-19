# Re-exports and wrappers for MOObject::multiOmicDataSet and related functions

#' `multiOmicDataSet` class
#'
#' Re-exported from [MOObject::multiOmicDataSet] -- view the MOObject docs for details.
#'
#' @inheritParams MOObject::multiOmicDataSet
#' @seealso [MOObject::multiOmicDataSet]
#' @importFrom MOObject multiOmicDataSet
#' @export
#' @family moo IO
multiOmicDataSet <- MOObject::multiOmicDataSet

#' Construct a multiOmicDataSet object from data frames
#'
#' Wraps [MOObject::create_multiOmicDataSet_from_dataframes()] to add a "colors" analysis if it is not already present.
#'
#' @inheritParams MOObject::create_multiOmicDataSet_from_dataframes
#' @seealso [MOObject::create_multiOmicDataSet_from_dataframes()]
#' @returns A [MOObject::multiOmicDataSet] object.
#' @export
#' @family moo IO
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
#' Wraps [MOObject::create_multiOmicDataSet_from_files()] to add a "colors" analysis if it is not already present.
#'
#' @inheritParams MOObject::create_multiOmicDataSet_from_files
#' @seealso [MOObject::create_multiOmicDataSet_from_files()]
#' @returns A [MOObject::multiOmicDataSet] object.
#' @export
#' @family moo IO
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
#' Re-exported from [MOObject::extract_counts] -- view the MOObject docs for details.
#'
#' @inheritParams MOObject::extract_counts
#' @seealso [MOObject::extract_counts()]
#' @returns A data frame of counts.
#' @importFrom MOObject extract_counts
#' @export
#' @family moo methods
extract_counts <- MOObject::extract_counts

#' Write a multiOmicDataSet to RDS
#'
#' Re-exported from [MOObject::write_multiOmicDataSet] -- view the MOObject docs for details.
#'
#' @inheritParams MOObject::write_multiOmicDataSet
#' @seealso [MOObject::write_multiOmicDataSet()]
#' @returns Invisibly returns `filepath`.
#' @importFrom MOObject write_multiOmicDataSet
#' @export
#' @family moo IO
write_multiOmicDataSet <- MOObject::write_multiOmicDataSet

#' Read a multiOmicDataSet from RDS
#'
#' Re-exported from [MOObject::read_multiOmicDataSet] -- view the MOObject docs for details.
#'
#' @inheritParams MOObject::read_multiOmicDataSet
#' @seealso [MOObject::read_multiOmicDataSet()]
#' @returns A [MOObject::multiOmicDataSet] object.
#' @importFrom MOObject read_multiOmicDataSet
#' @export
#' @family moo IO
read_multiOmicDataSet <- MOObject::read_multiOmicDataSet

#' Write multiOmicDataSet properties to individual files.
#'
#' Re-exported from [MOObject::write_multiOmicDataSet_properties] -- view the MOObject docs for details.
#'
#' @inheritParams MOObject::write_multiOmicDataSet_properties
#' @seealso [MOObject::write_multiOmicDataSet_properties()]
#' @returns Invisibly returns `output_dir`.
#' @importFrom MOObject write_multiOmicDataSet_properties
#' @export
#' @family moo IO
write_multiOmicDataSet_properties <- MOObject::write_multiOmicDataSet_properties

#' Read multiOmicDataSet properties from individual files.
#'
#' Re-exported from [MOObject::read_multiOmicDataSet_properties] -- view the MOObject docs for details.
#'
#' @inheritParams MOObject::read_multiOmicDataSet_properties
#' @seealso [MOObject::read_multiOmicDataSet_properties()]
#' @returns A [MOObject::multiOmicDataSet] object.
#' @importFrom MOObject read_multiOmicDataSet_properties
#' @export
#' @family moo IO
read_multiOmicDataSet_properties <- MOObject::read_multiOmicDataSet_properties
