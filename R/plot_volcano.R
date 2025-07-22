# plot_volcano_summary <- S7::new_generic("plot_volcano_summary", "moo_diff", function(moo_diff, ...) {
#   S7::S7_dispatch()
# })
#
# S7::method(plot_volcano_summary, multiOmicDataSet) <- function(moo_diff,
#                                                                ...) {
#   diff_dat <- moo_diff@analyses$diff %>% join_dfs_wide()
#   plot_volcano_summary(diff_dat, ...)
# }

# S7::method(plot_volcano_enhanced, multiOmicDataSet) <- function(moo_diff, ...) {
#   diff_dat <- moo_diff@analyses$diff %>% join_dfs_wide()
#   plot_volcano_enhanced(diff_dat, ...)
# }
