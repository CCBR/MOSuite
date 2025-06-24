# E2E workflow succeeds for RENEE data

    Code
      moo <- create_multiOmicDataSet_from_files(sample_meta_filepath = metadata_tsv,
        feature_counts_filepath = gene_counts_tsv) %>% clean_raw_counts() %>%
        filter_counts(group_colname = "condition", label_colname = "sample_id",
          minimum_count_value_to_be_considered_nonzero = 1,
          minimum_number_of_samples_with_nonzero_counts_in_total = 1,
          minimum_number_of_samples_with_nonzero_counts_in_a_group = 1, ) %>%
        normalize_counts(group_colname = "condition", label_colname = "sample_id") %>%
        diff_counts(covariates_colnames = "condition", contrast_colname = "condition",
          contrasts = c("knockout-wildtype"))
    Message
      Rows: 58929 Columns: 6
      -- Column specification --------------------------------------------------------
      Delimiter: "\t"
      chr (2): gene_id, GeneName
      dbl (4): KO_S3, KO_S4, WT_S1, WT_S2
      
      i Use `spec()` to retrieve the full column specification for this data.
      i Specify the column types or set `show_col_types = FALSE` to quiet this message.
      Rows: 4 Columns: 2
      -- Column specification --------------------------------------------------------
      Delimiter: "\t"
      chr (2): sample_id, condition
      
      i Use `spec()` to retrieve the full column specification for this data.
      i Specify the column types or set `show_col_types = FALSE` to quiet this message.
      * cleaning raw counts
      Not able to identify multiple id's in gene_id
      Columns that can be used to aggregate gene information gene_id
      Aggregating the counts for the same ID in different chromosome locations.
      Column used to Aggregate duplicate IDs: gene_id
      Number of rows before Collapse: 58929
      no duplicated IDs in gene_id
      * filtering clean counts
      Number of features after filtering: 291
      * normalizing filt counts
      Total number of features included: 291
      Sample columns: KO_S3, Sample columns: KO_S4, Sample columns: WT_S1, Sample columns: WT_S2
      * differential counts
      Setting first column of `counts` as gene annotation.
      Total number of genes included: 291

# E2E workflow succeeds for NIDAP data

    Code
      moo_nidap <- create_multiOmicDataSet_from_dataframes(sample_metadata = as.data.frame(
        nidap_sample_metadata), counts_dat = as.data.frame(nidap_raw_counts)) %>%
        clean_raw_counts() %>% filter_counts(group_colname = "Group") %>%
        normalize_counts(group_colname = "Group") %>% batch_correct_counts(
        covariates_colname = "Group", batch_colname = "Batch", label_colname = "Label") %>%
        diff_counts(count_type = "filt", sub_count_type = NULL, sample_id_colname = "Sample",
          feature_id_colname = "GeneName", covariates_colnames = c("Group", "Batch"),
          contrast_colname = c("Group"), contrasts = c("B-A", "C-A", "B-C"),
          input_in_log_counts = FALSE, return_mean_and_sd = TRUE,
          voom_normalization_method = "quantile", )
    Message
      * cleaning raw counts
      Not able to identify multiple id's in GeneName
      Columns that can be used to aggregate gene information GeneName
      Aggregating the counts for the same ID in different chromosome locations.
      Column used to Aggregate duplicate IDs: GeneName
      Number of rows before Collapse: 43280
      no duplicated IDs in GeneName
      * filtering clean counts
      Number of features after filtering: 7943
      * normalizing filt counts
      Total number of features included: 7943
      Sample columns: A1, Sample columns: A2, Sample columns: A3, Sample columns: B1, Sample columns: B2, Sample columns: B3, Sample columns: C1, Sample columns: C2, Sample columns: C3
      * batch-correcting norm-voom counts
      Found2batches
      Adjusting for2covariate(s) or covariate level(s)
      Standardizing Data across genes
      Fitting L/S model and finding priors
      Finding parametric adjustments
      Adjusting the Data
      
      The total number of features in output: 7943
      Number of samples after batch correction: 10
      * differential counts
      Setting first column of `counts` as gene annotation.
      Total number of genes included: 7942

