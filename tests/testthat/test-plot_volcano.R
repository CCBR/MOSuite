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
moo_nidap@analyses$diff <- nidap_deg_analysis_2
df_volc_sum <- plot_volcano_summary(nidap_deg_analysis)
df_volc_enh <- plot_volcano_enhanced(nidap_deg_analysis)
plot_volcano_enhanced(join_dfs(moo_nidap@analyses$diff))
