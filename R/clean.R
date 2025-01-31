# TODO replace print() with message() or warning() where appropriate

#' Clean Raw Counts
#'
#' @inheritParams filter_counts
#' @param cleanup_column_names Invalid raw counts column names can cause errors
#'   in the downstream analysis. If this is `TRUE`, any invalid column names
#'   will be automatically altered to a correct format. These format changes
#'   will include adding an "X" as the first character in any column name that
#'   began with a numeral and replacing some special characters ("-,:. ") with
#'   underscores ("_"). Invalid sample names and any changes made will be
#'   detailed in the template log.
#' @param split_gene_name If `TRUE`, split the gene name column by any of these special characters: `,|_-:`
#' @param aggregate_rows_with_duplicate_gene_names If a Feature ID (from the
#'   "Cleanup Column Names" parameter above) is found to be duplicated on
#'   multiple rows of the raw counts, the Log will report these Feature IDs.
#'   Using the default behavior (`TRUE`), the counts for all rows with a
#'   duplicate Feature IDs are aggregated into a single row. Counts are summed
#'   across duplicate Feature ID rows within each sample. Additional identifier
#'   columns, if present (e.g. Ensembl IDs), will be preserved and multiple
#'   matching identifiers in such additional columns will appear as
#'   comma-separated values in an aggregated row.
#' @param gene_name_column_to_use_for_collapsing_duplicates Select the column
#'   with Feature IDs to use as grouping elements to collapse the counts matrix.
#'   The log output will list the columns available to identify duplicate row
#'   IDs in order to aggregate information.
#'   If left blank your "Feature ID" Column will be used to Aggregate Rows. If
#'   "Feature ID" column can be split into multiple IDs the non Ensembl ID name
#'   will be used to aggregate duplicate IDs. If "Feature ID" column does not
#'   contain Ensembl IDs the split Feature IDs will be named 'Feature_id_1' and
#'   'Feature_id_2'. For this case an error will occur and you will have
#'   to manually enter the Column ID for this field.
#'
#' @returns `multiOmicDataSet` with cleaned counts
#' @export
#'
#' @examples
#' moo <- create_multiOmicDataSet_from_dataframes(
#'   as.data.frame(nidap_sample_metadata),
#'   as.data.frame(nidap_raw_counts),
#'   sample_id_colname = "Sample",
#' ) %>%
#'   clean_raw_counts(sample_id_colname = "Sample", feature_id_colname = "GeneName")
#' head(moo@counts$clean)
clean_raw_counts <- function(moo,
                             count_type = "raw",
                             sample_id_colname = NULL,
                             feature_id_colname = NULL,
                             samples_to_rename = c(""),
                             cleanup_column_names = TRUE,
                             split_gene_name = TRUE,
                             aggregate_rows_with_duplicate_gene_names = TRUE,
                             gene_name_column_to_use_for_collapsing_duplicates = "") {
  counts_dat <- moo@counts[[count_type]] %>% as.data.frame()
  sample_metadata <- moo@sample_meta %>% as.data.frame()

  if (is.null(sample_id_colname)) {
    sample_id_colname <- colnames(sample_metadata)[1]
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- colnames(counts_dat)[1]
  }
  # Sample Read Counts Plot
  read_plot <- plot_read_depth(counts_dat)

  # Manually rename samples
  counts_dat <- rename_samples(counts_dat, samples_to_rename)

  ### PH: START Clean up Sample Name columns
  ##################################
  ##### Cleanup Columns
  ##################################
  ## Look for any non standard Characters or automatic formatting introduced when Table is read included
  ## done after Rename step so that names do not automatically change before using names in Metadata.

  if (cleanup_column_names) {
    cl_og <- colnames(counts_dat)
    ## convert special charchers to _
    cl2 <- gsub("-| |\\:", "_", colnames(counts_dat))
    if (length(cl2[(cl2) != colnames(counts_dat)]) > 0) {
      print("Columns had special characters relpaced with _ ")
      # (colnames(counts_dat)[(colnames(counts_dat))!=cl2])
      # print(cl2[(cl2)!=colnames(counts_dat)])
      colnames(counts_dat) <- cl2
    }

    ## if names begin with number add X
    cl2 <- sub("^(\\d)", "X\\1", colnames(counts_dat))
    if (length(cl2[(cl2) != colnames(counts_dat)]) > 0) {
      print("Columns started with numbers and an X was added to colname :")
      # (colnames(counts_dat)[(colnames(counts_dat))!=cl2])
      # print(cl2[(cl2)!=colnames(counts_dat)])
      colnames(counts_dat) <- cl2
    }

    print(colnames(counts_dat)[!colnames(counts_dat) %in% feature_id_colname] %>% as.data.frame(),
      row.names = F
    )
    # print("Final Colnames:")
  } else {
    ## invalid name format
    if (any(make.names(colnames(counts_dat)) != colnames(counts_dat))) {
      print("Error: The following counts matrix column names are not valid:\n")
      print(colnames(counts_dat)[make.names(colnames(counts_dat)) != colnames(counts_dat)])
      print(
        "Likely causes are columns starting with numbers or other special characters eg spaces.\n"
      )
      # stop("Bad column names.")
    }
    ## Names Contain dashes
    if (sum(grepl("-", colnames(counts_dat))) != 0) {
      print("The sample names cannot contain dashes.")
      print(colnames(counts_dat)[grepl("-", colnames(counts_dat))])
      # stop("No dashes allowed in column names")
    }
  }
  ### PH: END Clean up Sample Name columns

  # Split Ensemble + Gene name
  counts_dat <- separate_gene_meta_columns(counts_dat,
    split_gene_name = split_gene_name
  )

  # Aggregate duplicate gene names
  counts_dat <- aggregate_duplicate_gene_names(counts_dat,
    gene_name_column_to_use_for_collapsing_duplicates = gene_name_column_to_use_for_collapsing_duplicates,
    aggregate_rows_with_duplicate_gene_names = aggregate_rows_with_duplicate_gene_names,
    split_gene_name = split_gene_name
  )

  moo@counts[["clean"]] <- counts_dat
  return(moo)
}

#' Remove version number from ENSEMBLE IDs
#'
#' @param x vector of IDs
#'
#' @return IDs without version numbers
#' @export
#'
#'
strip_ensembl_version <- function(x) {
  return(unlist(lapply(stringr::str_split(x, "[.]"), "[[", 1)))
}

#' Create read depth plot
#'
#' @param counts_dat dataframe with raw counts data
#'
#' @returns ggplot object
#' @export
#'
plot_read_depth <- function(counts_dat) {
  sample_names <- read_sums <- NULL

  sum_df <- counts_dat %>%
    dplyr::summarize(dplyr::across(tidyselect::where(is.numeric), sum)) %>%
    tidyr::pivot_longer(dplyr::everything(),
      names_to = "sample_names",
      values_to = "column_sums"
    )

  # Plotting
  read_plot <- ggplot2::ggplot(sum_df, ggplot2::aes(x = sample_names, y = read_sums)) +
    ggplot2::geom_bar(stat = "identity", fill = "blue") +
    ggplot2::labs(title = "Total Reads per Sample", x = "Samples", y = "Read Count") +
    ggplot2::scale_y_continuous(labels = scales::label_comma())
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
}

#' Separate gene metadata column
#'
#' @param counts_dat dataframe with raw counts data
#' @inheritParams clean_raw_counts
#'
#' @returns dataframe with metadata separated
separate_gene_meta_columns <- function(counts_dat, split_gene_name = TRUE) {
  ## Identify and separate Gene Name Columns into multiple Gene Metadata columns
  ##################################
  ## Split Ensemble + Gene name
  ##################################
  ## First check if Feature ID column  can be split by ",|_-:"
  ## Then check if one column contains Ensemble (regex '^ENS[A-Z]+[0-9]+')
  ##   check if Ensemble ID has version info and remove version
  ##   If one column contains Ensemble ID Assume other column is Gene names
  ## If Column does not contain Ensmeble ID name split columns Gene_ID_1 and Gene_ID_2

  ## if split_Gene_name ==F then will rename feature_id_colname column to either Gene(Bulk RNAseq) or FeatureID(Proteomics)

  feature_id_colname <- colnames(counts_dat)[1]
  out_colname <- feature_id_colname

  if (split_gene_name == T) {
    Ensembl_ID <- stringr::str_split_fixed(counts_dat[, feature_id_colname], "_|-|:|\\|", n = 2) %>%
      data.frame()
    EnsCol <- apply(Ensembl_ID, c(1, 2), function(x) {
      grepl("^ENS[A-Z]+[0-9]+", x)
    })

    if ("" %in% Ensembl_ID[, 1] | "" %in% Ensembl_ID[, 2]) {
      print(paste0("Not able to identify multiple id's in ", feature_id_colname))
      # colnames(df)[colnames(df)%in%clm]=gene_col
      colnames(counts_dat)[colnames(counts_dat) %in% feature_id_colname] <- out_colname
    } else {
      ## at least one column must have all ensemble ids found in EnsCol
      if (nrow(EnsCol[EnsCol[, 1] == T, ]) == nrow(Ensembl_ID) |
        nrow(EnsCol[EnsCol[, 2] == T, ]) == nrow(Ensembl_ID)) {
        colnames(Ensembl_ID)[colSums(EnsCol) != nrow(Ensembl_ID)] <- out_colname

        ## check if Ensmble column has version information
        if (grepl("^ENS[A-Z]+[0-9]+\\.[0-9]+$", Ensembl_ID[, colSums(EnsCol) == nrow(Ensembl_ID)]) %>% sum() == nrow(Ensembl_ID)) {
          colnames(Ensembl_ID)[colSums(EnsCol) == nrow(Ensembl_ID)] <- "Ensembl_ID_version"
          Ensembl_ID$Ensembl_ID <- strip_ensembl_version(Ensembl_ID$Ensembl_ID_version)
        } else {
          colnames(Ensembl_ID)[colSums(EnsCol) == nrow(Ensembl_ID)] <- "Ensembl_ID"
        }
      } else {
        colnames(Ensembl_ID) <- c("Feature_id_1", "Feature_id_2")
        print("Could not determine ID formats from split 'Feature ID' Column")
      }
      counts_dat <- cbind(Ensembl_ID, counts_dat[, !colnames(counts_dat) %in% feature_id_colname])
    }
  } else {
    colnames(counts_dat)[colnames(counts_dat) %in% feature_id_colname] <- out_colname
  }
  return(counts_dat)
}

#' Aggregate duplicate gene names
#'
#' @inheritParams clean_raw_counts
#' @inheritParams separate_gene_meta_columns
#'
#' @returns data frame with columns separated if possible
#'
aggregate_duplicate_gene_names <- function(counts_dat,
                                           gene_name_column_to_use_for_collapsing_duplicates,
                                           aggregate_rows_with_duplicate_gene_names,
                                           split_gene_name) {
  ##################################
  ## If duplicate gene, aggregate information to single row
  ##################################

  feature_id_colname <- colnames(counts_dat)[1]

  ## If user uses "Feature ID" column then switch to empty for appropriate behavior based on other parameters
  if (gene_name_column_to_use_for_collapsing_duplicates == feature_id_colname) {
    gene_name_column_to_use_for_collapsing_duplicates <- ""
  }
  if (gene_name_column_to_use_for_collapsing_duplicates == "" &
    ("Feature_id_1" %in% colnames(counts_dat)) == F) {
    gene_name_column_to_use_for_collapsing_duplicates <- feature_id_colname
  }

  ## Error Check if Column is Numeric
  nums <- unlist(lapply(counts_dat, is.numeric))
  nums <- names(nums[nums])
  print("")
  print("Columns that can be used to aggregate gene information")
  print(counts_dat[, !names(counts_dat) %in% nums, drop = F] %>% colnames())

  print("")

  ##########
  ## This section will Print duplicate row names when Aggregation column is not Specified.
  ## Purpose is to Identify Row Annotation columns and show user that rows may duplicated
  #######
  ## Print what rows are duplicated in selected annotation column
  ## if no column name given default to Gene(Bulk RNAseq) or FeatureID(Proteomics)
  ## Options:  1 Use default RowName for data_type
  ##           2 Use default Row name when Split Gene name recognizes Annotation name type
  ##           3 show duplicate rows for all row annotations columns when  Split Gene name does not recognize Annotation name type
  if (gene_name_column_to_use_for_collapsing_duplicates == "") {
    if (split_gene_name == F) {
      ## If no Column name given for Aggregation then display Feature ID duplicates
      print(paste0("genes with duplicate IDs in ", feature_id_colname, ":"))

      ## if Gene Name column is split then select Column Names generated from "Split Ensemble + Gene name"
      ## Raw If Feature_id_1 is generated it means that "Split Ensemble + Gene name" could not recognize Gene name format (EnsembleID or GeneName)
      ## and so default is to identify duplicicates in Feature_id_1 column
    } else if (split_gene_name == T &
      grepl("Feature_id_1", colnames(counts_dat)) == F) {
      print(paste0("genes with duplicate IDs in ", feature_id_colname, ":"))

      counts_dat[duplicated(counts_dat[, gene_name_column_to_use_for_collapsing_duplicates]), gene_name_column_to_use_for_collapsing_duplicates] %>%
        unique() %>%
        as.character() %>%
        write(stdout())
    } else if (split_gene_name == T &
      grepl("Feature_id_1", colnames(counts_dat)) == T) {
      print(paste0("genes with duplicate IDs in ", "Feature_id_1", ":"))

      counts_dat[duplicated(counts_dat[, "Feature_id_1"]), "Feature_id_1"] %>%
        unique() %>%
        as.character() %>%
        write(stdout())

      print(paste0("genes with duplicate IDs in ", "Feature_id_2", ":"))

      counts_dat[duplicated(counts_dat[, "Feature_id_2"]), "Feature_id_2"] %>%
        unique() %>%
        as.character() %>%
        write(stdout())
    }
  }

  ##########
  ## This section Aggregates duplicate Row names based on selected Annotation Column name
  #######
  if (aggregate_rows_with_duplicate_gene_names == TRUE) {
    print("Aggregating the counts for the same ID in different chromosome locations.")
    print("Column used to Aggregate duplicate IDs: ")
    print(gene_name_column_to_use_for_collapsing_duplicates)
    print("Number of rows before Collapse: ")
    print(nrow(counts_dat))

    if (sum(duplicated(counts_dat[, gene_name_column_to_use_for_collapsing_duplicates])) != 0) {
      print("")
      print("Duplicate IDs: ")
      print(counts_dat[duplicated(counts_dat[, gene_name_column_to_use_for_collapsing_duplicates]), gene_name_column_to_use_for_collapsing_duplicates] %>% as.character() %>% unique())

      dfagg <- counts_dat[, c(
        gene_name_column_to_use_for_collapsing_duplicates,
        nums
      )] %>%
        dplyr::group_by_at(gene_name_column_to_use_for_collapsing_duplicates) %>%
        dplyr::summarise_all(sum)

      if (ncol(counts_dat[, !names(counts_dat) %in% nums, drop = FALSE]) > 1) {
        ## collapse non-numeric columns
        dfagg2 <- counts_dat[, !names(counts_dat) %in% nums] %>%
          dplyr::group_by_at(gene_name_column_to_use_for_collapsing_duplicates) %>%
          dplyr::summarise_all(paste, collapse = ",")

        dfagg <- merge(
          dfagg2,
          dfagg,
          by = eval(gene_name_column_to_use_for_collapsing_duplicates),
          sort = FALSE
        ) %>% as.data.frame()
      }
      dfout <- dfagg
      print("Number of rows after Collapse: ")
      print(nrow(dfout))
    } else {
      print(
        paste0(
          "no duplicated IDs in ",
          gene_name_column_to_use_for_collapsing_duplicates
        )
      )
      dfout <- counts_dat
    }
  } else {
    if (gene_name_column_to_use_for_collapsing_duplicates != "") {
      print("")
      print(
        paste0(
          "Duplicate IDs in ",
          gene_name_column_to_use_for_collapsing_duplicates,
          " Column:"
        )
      )
      print(counts_dat[duplicated(counts_dat[, gene_name_column_to_use_for_collapsing_duplicates]), gene_name_column_to_use_for_collapsing_duplicates] %>% as.character() %>% unique())
    }

    print("")
    print(
      paste0(
        "If you desire to Aggregate row feature information select appropriate Column to use for collapsing duplicates"
      )
    )
  }

  return(counts_dat)
}
