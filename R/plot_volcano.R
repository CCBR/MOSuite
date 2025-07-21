plot_volcano_summary <- S7::new_generic("plot_volcano_summary", "moo_counts", function(moo_counts, ...) {
  S7::S7_dispatch()
})

S7::method(plot_volcano_summary, multiOmicDataSet) <- function(moo_counts,
                                                               count_type,
                                                               sub_count_type = NULL,
                                                               ...) {
  counts_dat <- extract_counts(moo_counts, count_type, sub_count_type)
  plot_volcano_summary(counts_dat, ...)
}
