test_that("E2E workflow succeeds for RENEE data", {
  options(moo_print_plots = FALSE, moo_save_plots = FALSE)
  gene_counts_tsv <- system.file(
    "extdata",
    "RSEM.genes.expected_count.all_samples.txt.gz",
    package = "MOSuite"
  )
  metadata_tsv <- system.file(
    "extdata",
    "sample_metadata.tsv.gz",
    package = "MOSuite"
  )

  expect_snapshot(
    moo <- create_multiOmicDataSet_from_files(
      sample_meta_filepath = metadata_tsv,
      feature_counts_filepath = gene_counts_tsv
    ) |>
      clean_raw_counts() |>
      filter_counts(
        group_colname = "condition",
        label_colname = "sample_id",
        minimum_count_value_to_be_considered_nonzero = 1,
        minimum_number_of_samples_with_nonzero_counts_in_total = 1,
        minimum_number_of_samples_with_nonzero_counts_in_a_group = 1,
      ) |>
      normalize_counts(
        group_colname = "condition",
        label_colname = "sample_id"
      ) |>
      diff_counts(
        covariates_colnames = "condition",
        contrast_colname = "condition",
        contrasts = c("knockout-wildtype")
      ) |>
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
  )
})

test_that("E2E workflow succeeds for NIDAP data", {
  options(moo_print_plots = FALSE, moo_save_plots = FALSE)
  expect_snapshot(
    moo_nidap <- create_multiOmicDataSet_from_dataframes(
      sample_metadata = as.data.frame(nidap_sample_metadata),
      counts_dat = as.data.frame(nidap_raw_counts)
    ) |>
      clean_raw_counts() |>
      filter_counts(group_colname = "Group") |>
      normalize_counts(group_colname = "Group") |>
      batch_correct_counts(
        covariates_colname = "Group",
        batch_colname = "Batch",
        label_colname = "Label"
      ) |>
      diff_counts(
        count_type = "filt",
        sub_count_type = NULL,
        sample_id_colname = "Sample",
        feature_id_colname = "GeneName",
        covariates_colnames = c("Group", "Batch"),
        contrast_colname = c("Group"),
        contrasts = c("B-A", "C-A", "B-C"),
        input_in_log_counts = FALSE,
        return_mean_and_sd = TRUE,
        voom_normalization_method = "quantile",
      ) |>
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
  )
})
