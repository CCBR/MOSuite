moo_nidap <- multiOmicDataSet(
  sample_metadata = as.data.frame(nidap_sample_metadata),
  anno_dat = data.frame(),
  counts_lst = list(
    "raw" = as.data.frame(nidap_raw_counts),
    "clean" = as.data.frame(nidap_clean_raw_counts),
    "filt" = as.data.frame(nidap_filtered_counts),
    "norm" = list("voom" = as.data.frame(nidap_norm_counts))
  )
)

test_that("differential analysis works for NIDAP", {
  options(moo_print_plots = FALSE)
  deg_moo <- moo_nidap %>%
    diff_counts(
      count_type = "filt",
      sub_count_type = NULL,
      sample_id_colname = "Sample",
      feature_id_colname = "Gene",
      covariates_colnames = c("Group", "Batch"),
      contrast_colname = c("Group"),
      contrasts = c("B-A", "C-A", "B-C"),
      voom_normalization_method = "quantile",
    )

  expect_equal(
    deg_moo@analyses$diff %>%
      join_dfs_wide() %>%
      dplyr::arrange(Gene) %>%
      dplyr::select(order(colnames(.))),
    nidap_deg_analysis_2 %>%
      join_dfs_wide() %>%
      dplyr::arrange(Gene) %>%
      dplyr::select(order(colnames(.))),
    tolerance = 0.01
  )
})

test_that("diff_counts works for RENEE on macOS", {
  skip_on_os("linux") # these expected values only work for macOS
  options(moo_print_plots = FALSE)
  moo_renee <- create_multiOmicDataSet_from_dataframes(
    readr::read_tsv(
      system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite")
    ),
    counts_dat = gene_counts
  ) %>%
    clean_raw_counts() %>%
    filter_counts(
      group_colname = "condition",
      label_colname = "sample_id",
      minimum_count_value_to_be_considered_nonzero = 1,
      minimum_number_of_samples_with_nonzero_counts_in_total = 1,
      minimum_number_of_samples_with_nonzero_counts_in_a_group = 1
    ) %>%
    normalize_counts(group_colname = "condition", label_colname = "sample_id")
  moo_renee <- moo_renee %>%
    diff_counts(
      count_type = "norm",
      sub_count_type = "voom",
      sample_id_colname = NULL,
      feature_id_colname = NULL,
      covariates_colnames = c("condition"),
      contrast_colname = c("condition"),
      # , 'condition2'), # TODO does not currently work for more than one contrast column
      contrasts = c("knockout-wildtype"),
      voom_normalization_method = "TMM",
      return_mean_and_sd = TRUE,
      input_in_log_counts = TRUE
    )
  actual <- moo_renee@analyses$diff[[1]] %>% head()
  expected <- structure(
    list(
      gene_id = c(
        "ENSG00000160179.18",
        "ENSG00000258017.1",
        "ENSG00000282393.1",
        "ENSG00000286104.1",
        "ENSG00000274422.1",
        "ENSG00000154734.15"
      ),
      knockout_mean = c(
        10.9805713961628,
        9.00423753343925,
        9.00423753343925,
        9.00423753343925,
        9.00423753343925,
        8.60833480887895
      ),
      knockout_sd = c(
        2.1542262015539,
        0.640731950886978,
        0.479050054020298,
        0.640731950886978,
        0.479050054020298,
        0.0808409484333391
      ),
      wildtype_mean = c(
        12.3499548758012,
        8.87501967069015,
        8.87501967069015,
        8.87501967069015,
        8.87501967069015,
        14.6328231986934
      ),
      wildtype_sd = c(
        0.082485020673847,
        0.00393703924789543,
        0.00393703924789543,
        0.00393703924789543,
        0.00393703924789543,
        0.00393703924789669
      ),
      FC = c(
        -2.54684822818721,
        1.11617198002536,
        1.07719571726849,
        1.11617198002536,
        1.07719571726849,
        -65.0581764060579
      ),
      logFC = c(
        -1.34871298907199,
        0.158559335064906,
        0.107280399074217,
        0.158559335064906,
        0.107280399074217,
        -6.02365847879717
      ),
      tstat = c(
        -1.27583132251236,
        0.432232276101188,
        0.344942879338297,
        0.432232276101188,
        0.344942879338297,
        -32.9652734167629
      ),
      pval = c(
        0.273465120057361,
        0.688646721521334,
        0.748130864876644,
        0.688646721521334,
        0.748130864876644,
        7.17809796209428e-06
      ),
      adjpval = c(
        0.506868470934344,
        0.79745817464873,
        0.79745817464873,
        0.79745817464873,
        0.79745817464873,
        0.00102459282397239
      )
    ),
    row.names = c(NA, 6L),
    class = "data.frame"
  )
  expect_equal(actual, expected)
})

test_that("diff_counts works for RENEE on linux", {
  skip_on_os("mac") # these expected values only work for linux
  options(moo_print_plots = FALSE)
  moo_renee <- create_multiOmicDataSet_from_dataframes(
    readr::read_tsv(
      system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite")
    ),
    counts_dat = gene_counts
  ) %>%
    clean_raw_counts() %>%
    filter_counts(
      group_colname = "condition",
      label_colname = "sample_id",
      minimum_count_value_to_be_considered_nonzero = 1,
      minimum_number_of_samples_with_nonzero_counts_in_total = 1,
      minimum_number_of_samples_with_nonzero_counts_in_a_group = 1
    ) %>%
    normalize_counts(group_colname = "condition", label_colname = "sample_id")
  moo_renee <- moo_renee %>%
    diff_counts(
      count_type = "norm",
      sub_count_type = "voom",
      sample_id_colname = NULL,
      feature_id_colname = NULL,
      covariates_colnames = c("condition"),
      contrast_colname = c("condition"),
      # , 'condition2'), # TODO does not currently work for more than one contrast column
      contrasts = c("knockout-wildtype"),
      voom_normalization_method = "TMM",
      return_mean_and_sd = TRUE,
      input_in_log_counts = TRUE
    )
  actual <- moo_renee@analyses$diff[[1]] %>% head()

  expected_linux <- structure(
    list(
      gene_id = c(
        "ENSG00000160179.18",
        "ENSG00000258017.1",
        "ENSG00000282393.1",
        "ENSG00000286104.1",
        "ENSG00000274422.1",
        "ENSG00000154734.15"
      ),
      knockout_mean = c(
        10.980894421743,
        9.00456055901947,
        9.00456055901947,
        9.00456055901947,
        9.00456055901947,
        8.60865783445917
      ),
      knockout_sd = c(
        2.15422620156956,
        0.640731950871317,
        0.479050054035958,
        0.640731950871317,
        0.479050054035958,
        0.0808409484176796
      ),
      wildtype_mean = c(
        12.3496318502785,
        8.87469664516751,
        8.87469664516751,
        8.87469664516751,
        8.87469664516751,
        14.6325001731707
      ),
      wildtype_sd = c(
        0.0815713665212029,
        0.00302338509525259,
        0.00302338509525259,
        0.00302338509525259,
        0.00302338509525259,
        0.00302338509525259
      ),
      FC = c(
        -2.5616561054721,
        1.16633707170805,
        1.04318046829749,
        1.16633707170805,
        1.04318046829749,
        -64.9400478688159
      ),
      logFC = c(
        -1.35707681126473,
        0.22198478802868,
        0.0609887630307995,
        0.22198478802868,
        0.0609887630307995,
        -6.02103654295574
      ),
      tstat = c(
        -1.28449015139078,
        0.69892909674399,
        0.225897873810806,
        0.69892909674399,
        0.225897873810806,
        -46.0999107324897
      ),
      pval = c(
        0.272572296938119,
        0.525504375492067,
        0.833052756338384,
        0.525504375492067,
        0.833052756338384,
        2.64449586643222e-06
      ),
      adjpval = c(
        0.486616800055169,
        0.698272937297678,
        0.877247422097368,
        0.698272937297678,
        0.877247422097368,
        0.000384774148565888
      )
    ),
    row.names = c(NA, 6L),
    class = "data.frame"
  )
  expect_equal(actual, expected_linux)
})

test_that("diff_counts errors", {
  expect_error(
    moo_nidap %>% diff_counts(count_type = "DoesNotExist"),
    "count_type DoesNotExist not in"
  )
  expect_error(
    moo_nidap %>%
      diff_counts(count_type = "raw", sub_count_type = "DoesNotExist"),
    "raw counts is not a named list"
  )
  expect_error(
    moo_nidap %>%
      diff_counts(count_type = "norm", sub_count_type = "DoesNotExist"),
    "sub_count_type DoesNotExist is not in"
  )
})

test_that("filter_diff works for NIDAP on macOS", {
  skip_on_os("linux")
  options(moo_print_plots = FALSE)
  moo <- moo_nidap %>%
    diff_counts(
      count_type = "filt",
      sub_count_type = NULL,
      sample_id_colname = "Sample",
      feature_id_colname = "Gene",
      covariates_colnames = c("Group", "Batch"),
      contrast_colname = c("Group"),
      contrasts = c("B-A", "C-A", "B-C"),
      voom_normalization_method = "quantile",
    ) %>%
    filter_diff(
      significance_column = "adjpval",
      significance_cutoff = 0.05,
      change_column = "logFC",
      change_cutoff = 1,
      filtering_mode = "any",
      include_estimates = c("FC", "logFC", "tstat", "pval", "adjpval"),
      round_estimates = TRUE,
      rounding_decimal_for_percent_cells = 0,
      contrast_filter = "none",
      contrasts = c(),
      groups = c(),
      groups_filter = "none",
      label_font_size = 6,
      label_distance = 1,
      y_axis_expansion = 0.08,
      fill_colors = c("steelblue1", "whitesmoke"),
      pie_chart_in_3d = TRUE,
      bar_width = 0.4,
      draw_bar_border = TRUE,
      plot_type = "bar",
      plot_titles_fontsize = 12
    )
  expect_equal(moo@analyses$diff_filt, nidap_deg_gene_list)
})

test_that("filter_diff works for NIDAP on linux", {
  skip_on_os("mac")
  options(moo_print_plots = FALSE)
  moo <- moo_nidap %>%
    diff_counts(
      count_type = "filt",
      sub_count_type = NULL,
      sample_id_colname = "Sample",
      feature_id_colname = "Gene",
      covariates_colnames = c("Group", "Batch"),
      contrast_colname = c("Group"),
      contrasts = c("B-A", "C-A", "B-C"),
      voom_normalization_method = "quantile",
    ) %>%
    filter_diff(
      significance_column = "adjpval",
      significance_cutoff = 0.05,
      change_column = "logFC",
      change_cutoff = 1,
      filtering_mode = "any",
      include_estimates = c("FC", "logFC", "tstat", "pval", "adjpval"),
      round_estimates = TRUE,
      rounding_decimal_for_percent_cells = 0,
      contrast_filter = "none",
      contrasts = c(),
      groups = c(),
      groups_filter = "none",
      label_font_size = 6,
      label_distance = 1,
      y_axis_expansion = 0.08,
      fill_colors = c("steelblue1", "whitesmoke"),
      pie_chart_in_3d = TRUE,
      bar_width = 0.4,
      draw_bar_border = TRUE,
      plot_type = "bar",
      plot_titles_fontsize = 12
    )
  expected_head <- structure(
    list(
      Gene = c(
        "1110034G24Rik",
        "3110082I17Rik",
        "4632428N05Rik",
        "4833439L19Rik",
        "4930523C07Rik",
        "5430427O19Rik"
      ),
      `B-A_FC` = c(21.7, -1.73, 2.43, -1.38, -2.3, -2.22),
      `B-A_logFC` = c(4.44, -0.789, 1.28, -0.46, -1.2, -1.15),
      `B-A_tstat` = c(3.2, -1.35, 2.76, -1.18, -1.62, -2.46),
      `B-A_pval` = c(0.00782, 0.203, 0.0177, 0.26, 0.133, 0.0307),
      `B-A_adjpval` = c(0.21, 0.71, 0.303, 0.758, 0.617, 0.377),
      `C-A_FC` = c(36.6, -21.9, 4.66, -3.59, 4.5, -4.49),
      `C-A_logFC` = c(5.2, -4.46, 2.22, -1.84, 2.17, -2.17),
      `C-A_tstat` = c(4.15, -3.8, 5.25, -3.76, 4.5, -3.68),
      `C-A_pval` = c(0.00141, 0.00265, 0.000222, 0.00281, 0.000767, 0.00327),
      `C-A_adjpval` = c(0.027, 0.0383, 0.00929, 0.0395, 0.0191, 0.0432),
      `B-C_FC` = c(-1.69, 12.7, -1.92, 2.61, -10.3, 2.02),
      `B-C_logFC` = c(-0.758, 3.67, -0.941, 1.38, -3.37, 1.01),
      `B-C_tstat` = c(-0.838, 2.93, -3.15, 2.63, -5.04, 1.53),
      `B-C_pval` = c(0.419, 0.0129, 0.00859, 0.0222, 0.000311, 0.153),
      `B-C_adjpval` = c(0.707, 0.144, 0.124, 0.186, 0.0224, 0.442)
    ),
    row.names = c(NA, 6L),
    class = "data.frame"
  )
  expected_tail <- structure(
    list(
      Gene = c("Zfand6", "Zfp35", "Zfp422", "Zfp706", "Zfp945", "Zhx1"),
      `B-A_FC` = c(1.22, 1.15, -1.43, -1.92, 10.5, -1.28),
      `B-A_logFC` = c(0.282, 0.198, -0.515, -0.938, 3.39, -0.356),
      `B-A_tstat` = c(0.808, 0.725, -1.92, -5.49, 2.94, -1.33),
      `B-A_pval` = c(0.435, 0.483, 0.0802, 0.00015, 0.0126, 0.209),
      `B-A_adjpval` = c(0.859, 0.871, 0.529, 0.0247, 0.258, 0.712),
      `C-A_FC` = c(2.19, -2.35, -2.39, -2.89, 21.5, 1.6),
      `C-A_logFC` = c(1.13, -1.23, -1.26, -1.53, 4.43, 0.68),
      `C-A_tstat` = c(3.61, -4.08, -4.61, -8.9, 4.26, 2.9),
      `C-A_pval` = c(0.00372, 0.00161, 0.000638, 1.47e-06, 0.00117, 0.0137),
      `C-A_adjpval` = c(0.0462, 0.0295, 0.0176, 0.000377, 0.0241, 0.0971),
      `B-C_FC` = c(-1.8, 2.69, 1.67, 1.51, -2.06, -2.05),
      `B-C_logFC` = c(-0.845, 1.43, 0.744, 0.594, -1.04, -1.04),
      `B-C_tstat` = c(-2.76, 4.59, 2.56, 3.16, -1.24, -4.23),
      `B-C_pval` = c(0.0175, 0.000656, 0.0253, 0.00845, 0.239, 0.00123),
      `B-C_adjpval` = c(0.165, 0.0329, 0.199, 0.123, 0.542, 0.0455)
    ),
    row.names = 630:635,
    class = "data.frame"
  )
  expect_equal(head(moo@analyses$diff_filt), expected_head)
  expect_equal(tail(moo@analyses$diff_filt), expected_tail)
})
