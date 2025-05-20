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
      dplyr::arrange(Gene) %>%
      dplyr::select(order(colnames(.))),
    nidap_deg_analysis_2 %>%
      dplyr::arrange(Gene) %>%
      dplyr::select(order(colnames(.)))
  )
})

test_that("diff_counts works for RENEE", {
  moo_renee <- create_multiOmicDataSet_from_dataframes(readr::read_tsv(
    system.file("extdata", "sample_metadata.tsv.gz", package = "MOSuite")
  ), gene_counts) %>%
    clean_raw_counts() %>%
    filter_counts(
      group_colname = "condition",
      label_colname = "sample_id",
      minimum_count_value_to_be_considered_nonzero = 1,
      minimum_number_of_samples_with_nonzero_counts_in_total = 1,
      minimum_number_of_samples_with_nonzero_counts_in_a_group = 1,
      print_plots = FALSE
    )
  moo_renee %<>%
    diff_counts(
      count_type = "filt",
      sub_count_type = NULL,
      sample_id_colname = NULL,
      feature_id_colname = NULL,
      covariates_colnames = c("condition"),
      contrast_colname = c("condition"),
      contrasts = c("knockout-wildtype"),
      voom_normalization_method = "TMM",
      return_mean_and_sd = TRUE,
    )
  expect_equal(
    head(moo_renee@analyses$diff),
    structure(
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
          11.7169579692258,
          9.98724215990714,
          9.98724215990714,
          9.98724215990714,
          9.98724215990714,
          9.19476090954656
        ),
        wildtype_mean = c(
          12.241401845636,
          8.43773668876116,
          8.43773668876116,
          8.43773668876116,
          8.43773668876116,
          15.3727745333725
        ),
        knockout_sd = c(
          1.48855413589981,
          0.957633420689776,
          1.28384204368286,
          0.957633420689776,
          1.28384204368286,
          0.163104311496542
        ),
        wildtype_sd = c(
          0.149946795889277,
          0.295929611328406,
          0.295929611328406,
          0.295929611328406,
          0.295929611328406,
          0.263927886911104
        ),
        `knockout-wildtype_FC` = c(
          -1.46050133433163,
          3.13355018539208,
          2.69559941120833,
          3.13355018539208,
          2.69559941120833, -72.4048093604422
        ),
        `knockout-wildtype_logFC` = c(
          -0.546463676231276,
          1.64779809880242,
          1.43060611580631,
          1.64779809880242,
          1.43060611580631, -6.17801362382593
        ),
        `knockout-wildtype_tstat` = c(
          -0.784757682830198,
          1.78109778854808,
          1.40983695907771,
          1.78109778854808,
          1.40983695907771, -7.82912473454161
        ),
        `knockout-wildtype_pval` = c(
          0.456059101505561,
          0.11422290979923,
          0.19767773950156,
          0.11422290979923,
          0.19767773950156,
          6.28737377569409e-05
        ),
        `knockout-wildtype_adjpval` = c(
          0.539899521348879,
          0.203919427923779,
          0.259118117995289,
          0.203919427923779,
          0.259118117995289,
          0.00203291752080776
        ),
        KO_S3 = c(
          10.6643912455677,
          10.6643912455677,
          9.07942874484659,
          10.6643912455677,
          9.07942874484659,
          9.07942874484659
        ),
        KO_S4 = c(
          12.7695246928838,
          9.31009307424653,
          10.8950555749677,
          9.31009307424653,
          10.8950555749677,
          9.31009307424653
        ),
        WT_S1 = c(
          12.3474302418265,
          8.64699052368537,
          8.64699052368537,
          8.64699052368537,
          8.64699052368537,
          15.1861493347934
        ),
        WT_S2 = c(
          12.1353734494455,
          8.22848285383694,
          8.22848285383694,
          8.22848285383694,
          8.22848285383694,
          15.5593997319516
        )
      ),
      row.names = c(NA, 6L),
      class = "data.frame"
    )
  )
  expect_equal(
    tail(moo_renee@analyses$diff),
    structure(
      list(
        gene_id = c(
          "ENSG00000157538.14",
          "ENSG00000160193.11",
          "ENSG00000182093.15",
          "ENSG00000182362.14",
          "ENSG00000173276.14",
          "ENSG00000237232.7"
        ),
        knockout_mean = c(
          12.5094392195864,
          10.7797234102677,
          11.3909196209359,
          10.3557249569902,
          12.328154179894,
          9.98724215990714
        ),
        wildtype_mean = c(
          12.2914162548016,
          10.3911819865654,
          11.4263766505111,
          8.43773668876116,
          11.7948594475942,
          8.43773668876116
        ),
        knockout_sd = c(
          0.367816403713495,
          0.163104311496542,
          0.701257658901503,
          1.47874678982587,
          0.624192165501767,
          1.28384204368286
        ),
        wildtype_sd = c(
          0.261621182826875,
          0.225183757807691,
          0.0395538495401341,
          0.295929611328406,
          0.481559519595963,
          0.295929611328406
        ),
        `knockout-wildtype_FC` = c(
          1.13637444669373,
          1.3152492302628, -1.01219564897455,
          3.79515730062804,
          1.38084877968178,
          2.69559941120833
        ),
        `knockout-wildtype_logFC` = c(
          0.184438295560613,
          0.395336205749295, -0.0174881779019369,
          1.92415968330423,
          0.465555334984355,
          1.43060611580631
        ),
        `knockout-wildtype_tstat` = c(
          0.384986798464922,
          0.639732354714475, -0.0283577389556018,
          1.74301758795661,
          0.823051165697135,
          1.40983695907771
        ),
        `knockout-wildtype_pval` = c(
          0.710662915437076,
          0.540908656919554,
          0.978097999873126,
          0.120962220563016,
          0.435226574685942,
          0.19767773950156
        ),
        `knockout-wildtype_adjpval` = c(
          0.774542728060633,
          0.617272232014079,
          0.981470751596827,
          0.204651198743242,
          0.523350963775244,
          0.259118117995289
        ),
        KO_S3 = c(
          12.2493537462889,
          10.6643912455677,
          11.8867836669042,
          11.401356839734,
          11.8867836669042,
          9.07942874484659
        ),
        KO_S4 = c(
          12.7695246928838,
          10.8950555749677,
          10.8950555749677,
          9.31009307424653,
          12.7695246928838,
          10.8950555749677
        ),
        WT_S1 = c(
          12.1064221423227,
          10.2319530244065,
          11.454345445743,
          8.64699052368537,
          11.454345445743,
          8.64699052368537
        ),
        WT_S2 = c(
          12.4764103672805,
          10.5504109487243,
          11.3984078552793,
          8.22848285383694,
          12.1353734494455,
          8.22848285383694
        )
      ),
      row.names = 285:290,
      class = "data.frame"
    )
  )
})

test_that("diff_counts errors", {
  expect_error(
    moo_nidap %>% diff_counts(count_type = "raw", sub_count_type = "DoesNotExist"),
    "raw counts is not a named list"
  )
  expect_error(
    moo_nidap %>% diff_counts(count_type = "norm", sub_count_type = "DoesNotExist"),
    "sub_count_type DoesNotExist is not in"
  )
})
