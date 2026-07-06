#' Plot read depth as a bar plot
#'
#' The first argument can be a `multiOmicDataset` object (`moo`) or a `data.frame` containing counts.
#' For a `moo`, choose which counts slot to use with `count_type` & (optionally) `sub_count_type`.
#'
#' @param moo_counts counts dataframe or `multiOmicDataSet` containing `count_type` & `sub_count_type` in the counts
#'   slot
#' @param ... arguments forwarded to method
#'
#' @return ggplot barplot
#'
#' @export
#' @examples
#' # multiOmicDataSet
#' moo <- multiOmicDataSet(
#'   sample_metadata = nidap_sample_metadata,
#'   anno_dat = data.frame(),
#'   counts_lst = list(
#'     "raw" = nidap_raw_counts,
#'     "clean" = nidap_clean_raw_counts
#'   )
#' )
#'
#' plot_read_depth(moo, count_type = "clean")
#'
#' # dataframe
#' plot_read_depth(nidap_clean_raw_counts)
#'
#' @details
#'
#' # Methods
#'
#' | link to docs  | class  |
#' |---|---|
#' | [plot_read_depth()] | `multiOmicDataSet` |
#' | [plot_read_depth()] | `data.frame`       |
#'
#' @seealso
#' - [plot_read_depth.multiOmicDataSet()]
#' - [plot_read_depth.data.frame()]
#'
#' @family plotters
#' @keywords plotters
#' @family moo methods
plot_read_depth <- S7::new_generic(
  "plot_read_depth",
  dispatch_args = "moo_counts"
)

#' Plot read depth for multiOmicDataSet
#'
#' @rdname plot_read_depth.multiOmicDataSet
#' @aliases plot_read_depth.multiOmicDataSet
#' @usage NULL
#'
#' @param count_type the type of counts to use. Must be a name in the counts slot (`names(moo@counts)`).
#' @param sub_count_type used if `count_type` is a list in the counts slot: specify the sub count type within the list.
#'   Must be a name in `names(moo@counts[[count_type]])`.
#' @param sample_id_colname column in sample metadata containing sample IDs.
#' @param group_colname sample metadata column used to color bars. Leave blank to use the current single-color bar fill.
#' @param color_values colors used when `group_colname` is supplied. Named vectors are matched to group values;
#'   unnamed vectors follow group order and are extended with MOSuite colors when too few colors are supplied.
#'
#' @return ggplot barplot
#'
#' @examples
#' # multiOmicDataSet
#' moo <- multiOmicDataSet(
#'   sample_metadata = nidap_sample_metadata,
#'   anno_dat = data.frame(),
#'   counts_lst = list(
#'     "raw" = nidap_raw_counts,
#'     "clean" = nidap_clean_raw_counts
#'   )
#' )
#'
#' plot_read_depth(moo, count_type = "clean")
#'
#' @seealso [plot_read_depth()] generic
#' @family plotters for multiOmicDataSets
S7::method(plot_read_depth, multiOmicDataSet) <- function(
  moo_counts,
  count_type,
  sub_count_type = NULL,
  sample_id_colname = NULL,
  group_colname = "",
  color_values = c(
    "#5954d6",
    "#e1562c",
    "#b80058",
    "#00c6f8",
    "#d163e6",
    "#00a76c",
    "#ff9287",
    "#008cf9",
    "#006e00",
    "#796880",
    "#FFA500",
    "#878500"
  ),
  ...
) {
  counts_dat <- extract_counts(moo_counts, count_type, sub_count_type)
  color_by_group <- !is.null(group_colname) && trimws(group_colname) != ""
  if (!isTRUE(color_by_group)) {
    return(plot_read_depth(counts_dat, ...))
  }

  return(plot_read_depth(
    counts_dat,
    sample_metadata = moo_counts@sample_meta,
    sample_id_colname = sample_id_colname,
    group_colname = group_colname,
    color_values = color_values,
    ...
  ))
}

#' Plot read depth for `data.frame`
#'
#' @rdname plot_read_depth.data.frame
#' @aliases plot_read_depth.data.frame
#' @usage NULL
#'
#' @param sample_metadata sample metadata dataframe, required when `group_colname` is supplied.
#' @param sample_id_colname column in sample metadata containing sample IDs.
#' @param group_colname sample metadata column used to color bars. Leave blank to use the current single-color bar fill.
#' @param color_values colors used when `group_colname` is supplied. Named vectors are matched to group values;
#'   unnamed vectors follow group order and are extended with MOSuite colors when too few colors are supplied.
#' @param ... additional arguments (ignored; accepted for compatibility with the moo dispatch)
#'
#' @return ggplot barplot
#'
#' @examples
#' # dataframe
#' plot_read_depth(nidap_clean_raw_counts)
#'
#' @seealso [plot_read_depth()] generic
#' @family plotters for counts dataframes
S7::method(plot_read_depth, S7::class_data.frame) <- function(
  moo_counts,
  sample_metadata = NULL,
  sample_id_colname = NULL,
  group_colname = "",
  color_values = c(
    "#5954d6",
    "#e1562c",
    "#b80058",
    "#00c6f8",
    "#d163e6",
    "#00a76c",
    "#ff9287",
    "#008cf9",
    "#006e00",
    "#796880",
    "#FFA500",
    "#878500"
  ),
  ...
) {
  sample_names <- column_sums <- NULL
  counts_dat <- moo_counts
  sum_df <- counts_dat |>
    dplyr::summarize(dplyr::across(tidyselect::where(is.numeric), sum)) |>
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to = "sample_names",
      values_to = "column_sums"
    )

  color_by_group <- !is.null(group_colname) && trimws(group_colname) != ""
  if (color_by_group) {
    if (is.null(sample_metadata)) {
      stop("sample_metadata is required when group_colname is supplied")
    }
    if (is.null(sample_id_colname)) {
      sample_id_colname <- colnames(sample_metadata)[1]
    }
    sum_df <- sum_df |>
      dplyr::left_join(
        sample_metadata,
        by = stats::setNames(sample_id_colname, "sample_names")
      ) |>
      dplyr::filter(!is.na(!!rlang::sym(group_colname)))
    color_values <- resolve_plot_colors(sum_df, group_colname, color_values)
    sum_df <- sum_df |>
      dplyr::mutate(
        !!rlang::sym(group_colname) := as.character(
          !!rlang::sym(group_colname)
        )
      )
  }

  # Plotting
  read_plot <- ggplot2::ggplot(
    sum_df,
    ggplot2::aes(x = sample_names, y = column_sums)
  )

  if (color_by_group) {
    read_plot <- read_plot +
      ggplot2::geom_bar(
        ggplot2::aes(fill = !!rlang::sym(group_colname)),
        stat = "identity"
      ) +
      ggplot2::scale_fill_manual(values = color_values)
  } else {
    read_plot <- read_plot +
      ggplot2::geom_bar(stat = "identity", fill = "blue")
  }

  read_plot <- read_plot +
    ggplot2::labs(
      title = "Total Reads per Sample",
      x = "Samples",
      y = "Read Count"
    ) +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 14
      ),
      axis.text.y = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 16),
      plot.title = ggplot2::element_text(size = 20)
    )
  return(read_plot)
}
