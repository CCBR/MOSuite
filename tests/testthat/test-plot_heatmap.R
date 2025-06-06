colors_vec <- c(
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
)
test_that("correlation heatmap works", {
  p <- plot_corr_heatmap(
    nidap_filtered_counts %>%
      dplyr::select(tidyselect::all_of(
        c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")
      )) %>%
      as.data.frame(),
    sample_metadata = as.data.frame(nidap_sample_metadata),
    sample_id_colname = "Sample",
    label_colname = "Label",
    group_colname = "Group",
    color_values = colors_vec
  )
  expect_s4_class(p, "Heatmap")
  expect_equal(p@matrix, structure(
    c(
      0,
      0.0349436686640265,
      0.0346707459478316,
      0.136436292332774,
      0.145504435322238,
      0.165715414451367,
      0.266892740192968,
      0.310669867608233,
      0.281044606820525,
      0.0349436686640265,
      0,
      0.026894587617103,
      0.139740412452416,
      0.144024768162842,
      0.156995761981556,
      0.256516858290949,
      0.289312631501573,
      0.275306551391499,
      0.0346707459478316,
      0.026894587617103,
      0,
      0.113462670330904,
      0.132707949438905,
      0.146137196349944,
      0.2518982836951,
      0.307893394909586,
      0.277982555134354,
      0.136436292332774,
      0.139740412452416,
      0.113462670330904,
      0,
      0.0467104077868874,
      0.0778256905442303,
      0.18935488124329,
      0.238494284141649,
      0.209007629325352,
      0.145504435322238,
      0.144024768162842,
      0.132707949438905,
      0.0467104077868874,
      0,
      0.0532124359450156,
      0.140242067145314,
      0.179723372754429,
      0.15251602311055,
      0.165715414451367,
      0.156995761981556,
      0.146137196349944,
      0.0778256905442303,
      0.0532124359450156,
      0,
      0.141067943113981,
      0.160160263560895,
      0.14605974755951,
      0.266892740192968,
      0.256516858290949,
      0.2518982836951,
      0.18935488124329,
      0.140242067145314,
      0.141067943113981,
      0,
      0.104501003317621,
      0.0500950924722408,
      0.310669867608233,
      0.289312631501573,
      0.307893394909586,
      0.238494284141649,
      0.179723372754429,
      0.160160263560895,
      0.104501003317621,
      0,
      0.0899444885709063,
      0.281044606820525,
      0.275306551391499,
      0.277982555134354,
      0.209007629325352,
      0.15251602311055,
      0.14605974755951,
      0.0500950924722408,
      0.0899444885709063,
      0
    ),
    dim = c(9L, 9L),
    dimnames = list(
      c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"),
      c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")
    )
  ))
})

test_that("plot_corr_heatmap method dispatch works", {
  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts),
      "norm" = list(
        "voom" = as.data.frame(nidap_norm_counts)
      )
    )
  )
  expect_equal(
    plot_corr_heatmap(moo, "filt")@matrix,
    plot_corr_heatmap(
      as.data.frame(nidap_filtered_counts),
      sample_metadata = as.data.frame(nidap_sample_metadata),
      feature_id_colname = "Gene"
    )@matrix
  )
})

# TODO get heatmap working on tibbles also
# test_that("heatmap works", {
#   corHM <- plot_corr_heatmap(
#     counts_dat = nidap_filtered_counts %>%
#       dplyr::select(tidyselect::all_of(c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"))),
#     sample_metadata = nidap_sample_metadata,
#     sample_id_colname = "Sample",
#     label_colname = "Label",
#     group_column = "Group",
#     color_values = colors_vec
#   )
# })

test_that("plot_expr_heatmap works", {
  moo <- multiOmicDataSet(
    sample_metadata = as.data.frame(nidap_sample_metadata),
    anno_dat = data.frame(),
    counts_lst = list(
      "raw" = as.data.frame(nidap_raw_counts),
      "clean" = as.data.frame(nidap_clean_raw_counts),
      "filt" = as.data.frame(nidap_filtered_counts),
      "norm" = list(
        "voom" = as.data.frame(nidap_norm_counts)
      )
    )
  )
  # TODO refactor plot_expr_heatmap to avoid triggering warnings
  expect_warning(
    expect_warning(
      expect_warning({
        set.seed(20250226)
        p_moo <- plot_expr_heatmap(moo,
          count_type = "norm",
          sub_count_type = "voom",
          feature_id_colname = "Gene"
        )
      })
    )
  )
  expect_equal(
    head(p_moo@matrix),
    structure(c(
      -0.469646298150098, -0.999122775178692, -1.63914985230916,
      -1.11065487632401, 1.06508157716121, 0.129357631742822, -1.6308124461643,
      -0.700059853143015, -1.5189464375787, -0.224167735591629, 1.06508157716121,
      -1.56168918461362, -1.6308124461643, -1.89235880186714, -0.541702967261003,
      -2.0463972593436, 1.06508157716121, -1.79694425364586, 0.329750044812768,
      0.0301570520160167, 0.167894484181845, -0.0423239634814564, -0.120254965774652,
      0.199342022111398, 0.642860459417498, 0.367952795408648, 0.373030023826938,
      0.474515861027971, -0.0390725679377215, 0.117730807477454, 0.728033358831939,
      0.385380968911897, 0.599797184786654, 0.453174107982841, -0.104093868102577,
      0.461536877004385, 0.623852071863989, 0.835927627565177, 0.888284230577201,
      0.725967683142828, -1.99037275585028, 0.831065997888202, 0.825036578059704,
      1.0208559858833, 0.853098375022516, 0.922704248230578, -0.758155232754779,
      0.657242141825761, 0.581738677492799, 0.951267000403811, 0.817694958753715,
      0.847181934356472, -0.183295341063619, 0.962357960209459
    ), dim = c(
      6L,
      9L
    ), dimnames = list(c(
      "Il2rb", "Rora", "Tcrg-C1", "Pdcd1", "Dntt",
      "Eya2"
    ), c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3")))
  )

  expect_warning(expect_warning(expect_warning({
    set.seed(20250226)
    p_dat <- plot_expr_heatmap(
      as.data.frame(nidap_norm_counts),
      sample_metadata = as.data.frame(nidap_sample_metadata),
      feature_id_colname = "Gene"
    )
  })))

  expect_equal(
    p_moo@matrix,
    p_dat@matrix
  )
})
