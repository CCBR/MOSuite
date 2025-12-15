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

  # x <- nidap_deg_analysis %>%
  #   dplyr::arrange(Gene) %>%
  #   dplyr::select(order(colnames(.))) %>% dplyr::select(-C2) %>%
  #   as.data.frame()
  # y <- moo@analyses$diff %>%
  #   dplyr::arrange(Gene) %>%
  #   dplyr::select(order(colnames(.)))
  # equal_dfs(x, y)

  expect_equal(
    deg_moo@analyses$diff %>%
      join_dfs_wide() %>%
      dplyr::arrange(Gene) %>%
      dplyr::select(order(colnames(.))),
    nidap_deg_analysis_2 %>%
      join_dfs_wide() %>%
      dplyr::arrange(Gene) %>%
      dplyr::select(order(colnames(.)))
  )
})

test_that("diff_counts works for RENEE", {
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
  actual <- moo_renee@analyses$diff[[1]] %>%
    head() %>%
    dplyr::mutate(dplyr::across(
      dplyr::where(is.numeric),
      ~ round(.x, digits = 0)
    ))
  expected <- tibble::tibble(
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
      0.479050054020299,
      0.640731950886978,
      0.479050054020299,
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
    `FC` = c(
      -2.56266458541781,
      1.16581746593108,
      1.04271170723614,
      1.16581746593108,
      1.04271170723614,
      -64.968995694703
    ),
    `logFC` = c(
      -1.35764466371701,
      0.221341920947646,
      0.0603403313641113,
      0.221341920947646,
      0.0603403313641113,
      -6.02167949874076
    ),
    `tstat` = c(
      -1.28601115515934,
      0.696457513565273,
      0.223206820569603,
      0.696457513565273,
      0.223206820569603,
      -45.9646621635613
    ),
    `pval` = c(
      0.271954308769351,
      0.526807064335256,
      0.83497735832889,
      0.526807064335256,
      0.83497735832889,
      2.61647311409974e-06
    ),
    `adjpval` = c(
      0.485513520563689,
      0.700003907404382,
      0.878543258151392,
      0.700003907404382,
      0.878543258151392,
      0.000380696838101513
    )
  ) %>%
    dplyr::mutate(dplyr::across(
      dplyr::where(is.numeric),
      ~ round(.x, digits = 0)
    )) %>%
    as.data.frame()
  expect_equal(actual, expected)
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

test_that("filter_diff works for NIDAP", {
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
