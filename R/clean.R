# Clean Raw Counts [CCBR] (5453b016-53cf-44c8-b09b-7efda66543af): v82
#' @export
clean_raw_counts <- function(moo,
                             sample_names_column = "Sample",
                             gene_names_column = "GeneName",
                             samples_to_rename = c(""),
                             cleanup_column_names = TRUE,
                             split_gene_name = TRUE,
                             aggregate_rows_with_duplicate_gene_names = TRUE,
                             gene_name_column_to_use_for_collapsing_duplicates = "",
                             data_type = "Bulk RNAseq" # TODO refactor so this param isn't needed
) {
  # TODO delete library statements, use pkg::fcn syntax
  library(stringr)
  library(tidyr)
  library(dplyr)
  library(ggplot2)

  raw_counts_matrix <- moo@counts[["raw"]]
  sample_metadata <- moo@sample_meta

  # Sample Read Counts Plot
  read_plot <- plot_read_depth(raw_counts_matrix)

  ### PH: START Sample Validation
  ##### Sample Name Validation
  ## check if sample names are different between raw counts
  ## and metadata tables
  # TODO move this to S7 validator
  check_sample_names(raw_counts_matrix, sample_metadata, sample_names_column)
  ##### Sample Name Check for duplicated
  ## duplicate col name
  # TODO move this to S7 validator
  if (sum(duplicated(colnames(raw_counts_matrix))) != 0) {
    print("Duplicate column names are not allowed, the following columns were duplicated.\n")
    colnames(raw_counts_matrix)[duplicated(colnames(raw_counts_matrix))]
    stop("Duplicated columns")
  }
  ### PH: END Sample Validation


  # Manually rename samples
  raw_counts_matrix <- rename_samples(raw_counts_matrix, samples_to_rename)



  ### PH: START Clean up Sample Name columns
  ##################################
  ##### Cleanup Columns
  ##################################
  ## Look for any non standard Characters or automatic formatting introduced when Table is read included
  ## done after Rename step so that names do not automatically change before using names in Metadata.

  if (cleanup_column_names) {
    cl_og <- colnames(raw_counts_matrix)
    ## convert special charchers to _
    cl2 <- gsub("-| |\\:", "_", colnames(raw_counts_matrix))
    if (length(cl2[(cl2) != colnames(raw_counts_matrix)]) > 0) {
      print("Columns had special characters relpaced with _ ")
      # (colnames(raw_counts_matrix)[(colnames(raw_counts_matrix))!=cl2])
      # print(cl2[(cl2)!=colnames(raw_counts_matrix)])
      colnames(raw_counts_matrix) <- cl2
    }

    ## if names begin with number add X
    cl2 <- sub("^(\\d)", "X\\1", colnames(raw_counts_matrix))
    if (length(cl2[(cl2) != colnames(raw_counts_matrix)]) > 0) {
      print("Columns started with numbers and an X was added to colname :")
      # (colnames(raw_counts_matrix)[(colnames(raw_counts_matrix))!=cl2])
      # print(cl2[(cl2)!=colnames(raw_counts_matrix)])
      colnames(raw_counts_matrix) <- cl2
    }

    print(colnames(raw_counts_matrix)[!colnames(raw_counts_matrix) %in% gene_names_column] %>% as.data.frame(),
      row.names = F
    )
    # print("Final Colnames:")
  } else {
    ## invalid name format
    if (any(make.names(colnames(raw_counts_matrix)) != colnames(raw_counts_matrix))) {
      print("Error: The following counts matrix column names are not valid:\n")
      print(colnames(raw_counts_matrix)[make.names(colnames(raw_counts_matrix)) != colnames(raw_counts_matrix)])
      print(
        "Likely causes are columns starting with numbers or other special characters eg spaces.\n"
      )
      # stop("Bad column names.")
    }
    ## Names Contain dashes
    if (sum(grepl("-", colnames(raw_counts_matrix))) != 0) {
      print("The sample names cannot contain dashes.")
      print(colnames(raw_counts_matrix)[grepl("-", colnames(raw_counts_matrix))])
      # stop("No dashes allowed in column names")
    }
  }
  ### PH: END Clean up Sample Name columns

  # Split Ensemble + Gene name
  raw_counts_matrix <- separate_gene_meta_columns(raw_counts_matrix,
    gene_names_column = gene_names_column,
    split_gene_name = split_gene_name,
    data_type = data_type
  )

  # Aggregate duplicate gene names
  raw_counts_matrix <- aggregate_duplicate_gene_names(raw_counts_matrix,
    gene_names_column = gene_names_column,
    gene_name_column_to_use_for_collapsing_duplicates = gene_name_column_to_use_for_collapsing_duplicates,
    data_type = data_type,
    aggregate_rows_with_duplicate_gene_names = aggregate_rows_with_duplicate_gene_names
  )

  print(data_type)
  print(read_plot)

  moo@counts[["clean"]] <- raw_counts_matrix
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

plot_read_depth <- function(raw_counts_matrix) {
  # Exclude the gene column from sample columns for counts
  # TODO: do not assume the first column is the gene column
  sample_cols <- raw_counts_matrix[, -1]

  # Sum of each sample column
  column_sums <- colSums(sample_cols, na.rm = TRUE)
  # Create a data frame for plotting
  sum_df <- data.frame(
    sample_names = names(column_sums),
    read_sums = column_sums
  )

  # Plotting
  read_plot <- ggplot(sum_df, aes(x = sample_names, y = read_sums)) +
    geom_bar(stat = "identity", fill = "blue") +
    labs(title = "Total Reads per Sample", x = "Samples", y = "Read Count") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 14
      ),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 20)
    )
}

check_sample_names <- function(counts, metadata, sample_names_column) {
  raw_count_names <- colnames(counts)
  metadata_names <- metadata[, sample_names_column]

  different_names <- setdiff(metadata_names, raw_count_names)
  if (length(different_names) > 0) {
    stop(
      "The following sample names are different in the metadata but not the raw counts: ",
      paste(different_names, collapse = ",")
    )
  }
}

separate_gene_meta_columns <- function(raw_counts_matrix, gene_names_column = "GeneName", split_gene_name = TRUE, data_type = "Bulk RNAseq") {
  ## Identify and separate Gene Name Columns into multiple Gene Metadata columns
  ##################################
  ## Split Ensemble + Gene name
  ##################################
  ## First check if Feature ID column  can be split by ",|_-:"
  ## Then check if one column contains Ensemble (regex '^ENS[A-Z]+[0-9]+')
  ##   check if Ensemble ID has version info and remove version
  ##   If one column contains Ensemble ID Assume other column is Gene names
  ## If Column does not contain Ensmeble ID name split columns Gene_ID_1 and Gene_ID_2

  ## if split_Gene_name ==F then will rename gene_names_column column to either Gene(Bulk RNAseq) or FeatureID(Proteomics)

  print("")

  if (split_gene_name == T) {
    Ensembl_ID <- str_split_fixed(raw_counts_matrix[, gene_names_column], "_|-|:|\\|", n = 2) %>% data.frame()
    EnsCol <- apply(Ensembl_ID, c(1, 2), function(x) {
      grepl("^ENS[A-Z]+[0-9]+", x)
    })

    if ("" %in% Ensembl_ID[, 1] | "" %in% Ensembl_ID[, 2]) {
      print(paste0("Not able to identify multiple id's in ", gene_names_column))
      # colnames(df)[colnames(df)%in%clm]=gene_col
      if (data_type == "Bulk RNAseq") {
        colnames(raw_counts_matrix)[colnames(raw_counts_matrix) %in% gene_names_column] <- "Gene"
      } else if (data_type == "Proteomics") {
        colnames(raw_counts_matrix)[colnames(raw_counts_matrix) %in% gene_names_column] <- "FeatureID"
      } else {
        print("incorrect Data Type")
        incorrect_data_type
      }
    } else {
      ## at least one column must have all ensemble ids found in EnsCol
      if (nrow(EnsCol[EnsCol[, 1] == T, ]) == nrow(Ensembl_ID) |
        nrow(EnsCol[EnsCol[, 2] == T, ]) == nrow(Ensembl_ID)) {
        if (data_type == "Bulk RNAseq") {
          colnames(Ensembl_ID)[colSums(EnsCol) != nrow(Ensembl_ID)] <- "Gene"
        } else if (data_type == "Proteomics") {
          colnames(Ensembl_ID)[colSums(EnsCol) != nrow(Ensembl_ID)] <- "FeatureID"
        }
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
      raw_counts_matrix <- cbind(Ensembl_ID, raw_counts_matrix[, !colnames(raw_counts_matrix) %in% gene_names_column])
    }
  } else {
    if (data_type == "Bulk RNAseq") {
      colnames(raw_counts_matrix)[colnames(raw_counts_matrix) %in% gene_names_column] <- "Gene"
    } else if (data_type == "Proteomics") {
      colnames(raw_counts_matrix)[colnames(raw_counts_matrix) %in% gene_names_column] <- "FeatureID"
    } else {
      print("incorrect Data Type")
      incorrect_data_type
    }
  }
  return(raw_counts_matrix)
}

aggregate_duplicate_gene_names <- function(raw_counts_matrix, gene_names_column,
                                           gene_name_column_to_use_for_collapsing_duplicates,
                                           data_type,
                                           aggregate_rows_with_duplicate_gene_names) {
  ##################################
  ## If duplicate gene, aggregate information to single row
  ##################################

  ## If user uses "Feature ID" column then switch to empty for appropriate behavior based on other parameters
  if (gene_name_column_to_use_for_collapsing_duplicates == gene_names_column) {
    gene_name_column_to_use_for_collapsing_duplicates <- ""
  }
  ## Use different Gene Names Column Name based on data Type
  ## I think we should deprecate this and just stick with "FeatureID"
  if (gene_name_column_to_use_for_collapsing_duplicates == "" &
    ("Feature_id_1" %in% colnames(raw_counts_matrix)) == F) {
    if (data_type == "Bulk RNAseq") {
      gene_name_column_to_use_for_collapsing_duplicates <- "Gene"
    } else if (data_type == "Proteomics") {
      gene_name_column_to_use_for_collapsing_duplicates <- "FeatureID"
    }
  }

  ## Error Check if Column is Numeric
  nums <- unlist(lapply(raw_counts_matrix, is.numeric))
  nums <- names(nums[nums])
  print("")
  print("Columns that can be used to aggregate gene information")
  print(raw_counts_matrix[, !names(raw_counts_matrix) %in% nums, drop = F] %>% colnames())

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
      print(paste0("genes with duplicate IDs in ", gene_names_column, ":"))

      ## Print original Column name for user Reference then use new Column name to subset table
      if (data_type == "Bulk RNAseq") {
        gene_names_column <- "Gene"
      } else if (data_type == "Proteomics") {
        gene_names_column <- "FeatureID"
      }
      raw_counts_matrix[duplicated(raw_counts_matrix[, gene_names_column]), gene_names_column] %>%
        unique() %>%
        as.character() %>%
        write(stdout())

      ## if Gene Name column is split then select Column Names generated from "Split Ensemble + Gene name"
      ## Raw If Feature_id_1 is generated it means that "Split Ensemble + Gene name" could not recognize Gene name format (EnsembleID or GeneName)
      ## and so default is to identify duplicicates in Feature_id_1 column
    } else if (split_gene_name == T &
      grepl("Feature_id_1", colnames(raw_counts_matrix)) == F) {
      if (data_type == "Bulk RNAseq") {
        gene_names_column <- "Gene"
      } else if (data_type == "Proteomics") {
        gene_names_column <- "FeatureID"
      }
      print(paste0("genes with duplicate IDs in ", gene_names_column, ":"))

      raw_counts_matrix[duplicated(raw_counts_matrix[, gene_name_column_to_use_for_collapsing_duplicates]), gene_name_column_to_use_for_collapsing_duplicates] %>%
        unique() %>%
        as.character() %>%
        write(stdout())
    } else if (split_gene_name == T &
      grepl("Feature_id_1", colnames(raw_counts_matrix)) == T) {
      print(paste0("genes with duplicate IDs in ", "Feature_id_1", ":"))

      raw_counts_matrix[duplicated(raw_counts_matrix[, "Feature_id_1"]), "Feature_id_1"] %>%
        unique() %>%
        as.character() %>%
        write(stdout())

      print(paste0("genes with duplicate IDs in ", "Feature_id_2", ":"))

      raw_counts_matrix[duplicated(raw_counts_matrix[, "Feature_id_2"]), "Feature_id_2"] %>%
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
    print(nrow(raw_counts_matrix))

    if (sum(duplicated(raw_counts_matrix[, gene_name_column_to_use_for_collapsing_duplicates])) != 0) {
      print("")
      print("Duplicate IDs: ")
      print(raw_counts_matrix[duplicated(raw_counts_matrix[, gene_name_column_to_use_for_collapsing_duplicates]), gene_name_column_to_use_for_collapsing_duplicates] %>% as.character() %>% unique())

      dfagg <- raw_counts_matrix[, c(
        gene_name_column_to_use_for_collapsing_duplicates,
        nums
      )] %>%
        group_by_at(gene_name_column_to_use_for_collapsing_duplicates) %>%
        summarise_all(sum)

      if (ncol(raw_counts_matrix[, !names(raw_counts_matrix) %in% nums, drop = FALSE]) > 1) {
        ## collapse non-numeric columns
        dfagg2 <- raw_counts_matrix[, !names(raw_counts_matrix) %in% nums] %>%
          group_by_at(gene_name_column_to_use_for_collapsing_duplicates) %>%
          summarise_all(paste, collapse = ",")

        dfagg <- merge(
          dfagg2,
          dfagg,
          by = eval(gene_name_column_to_use_for_collapsing_duplicates),
          sort = F
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
      dfout <- raw_counts_matrix
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
      print(raw_counts_matrix[duplicated(raw_counts_matrix[, gene_name_column_to_use_for_collapsing_duplicates]), gene_name_column_to_use_for_collapsing_duplicates] %>% as.character() %>% unique())
    }

    print("")
    print(
      paste0(
        "If you desire to Aggregate row feature information select appropriate Column to use for collapsing duplicates"
      )
    )
  }

  return(raw_counts_matrix)
}
