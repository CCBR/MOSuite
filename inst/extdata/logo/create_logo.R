#!/usr/bin/env Rscript
# Create hex logo for MOSuite package - Version 4
# Design: Combine volcano plot + heatmap

library(ggplot2)
library(hexSticker)
library(ggimage)
library(magick)
library(tidyverse)

# Create directory if it doesn't exist
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)

set.seed(123)

# Emoji watermark layer - B&W cow emoji inside O of MOSuite
emoji_url <- "https://raw.githubusercontent.com/twitter/twemoji/master/assets/72x72/1f42e.png"
emoji_temp <- tempfile(fileext = ".png")
download.file(emoji_url, emoji_temp, quiet = TRUE, mode = "wb")

# Convert to grayscale and apply transparency
emoji_gray <- magick::image_read(emoji_temp)
emoji_gray <- magick::image_quantize(emoji_gray, colorspace = "gray")
# Apply transparency via channel operations
#emoji_gray <- magick::image_transparent(emoji_gray, color = "white", fuzz = 30)
emoji_gray <- magick::image_fx(emoji_gray, "0.7*u", channel = "alpha")
magick::image_write(emoji_gray, emoji_temp)

emoji_df <- data.frame(
  x = 0.366,
  y = 0.73,
  image = emoji_temp
)

emoji_layer <- ggimage::geom_image(
  data = emoji_df,
  aes(x = x, y = y, image = image),
  size = 0.06,
  inherit.aes = FALSE
)

# ---- Volcano plot layer ----
# Continuous V-shape distribution
# n_genes <- 600
# log2fc <- rnorm(n_genes, mean = 0, sd = 1.6)
# neg_log10_pval <- abs(log2fc) * 1.2 + rnorm(n_genes, mean = 0.4, sd = 0.4)
# neg_log10_pval <- pmax(neg_log10_pval, 0)
# volcano_data <- data.frame(
#   log2fc = log2fc,
#   neg_log10_pval = neg_log10_pval
# )

load(here::here('data', 'nidap_deg_analysis.rda'))
volcano_data <- nidap_deg_analysis |>
  rename(log2fc = `C-A_logFC`) |>
  mutate(neg_log10_pval = -log10(`C-A_pval`))

fc_cut <- 1.0
p_cut <- 2.5

volcano_data$regulation <- "Not significant"
volcano_data$regulation[volcano_data$neg_log10_pval >= p_cut] <- "Significant"
volcano_data$regulation[
  volcano_data$log2fc >= fc_cut & volcano_data$neg_log10_pval >= p_cut
] <- "Upregulated"
volcano_data$regulation[
  volcano_data$log2fc <= -fc_cut & volcano_data$neg_log10_pval >= p_cut
] <- "Downregulated"

volcano_data$x <- scales::rescale(volcano_data$log2fc, to = c(0.08, 0.92))
volcano_data$y <- scales::rescale(
  volcano_data$neg_log10_pval,
  to = c(0.10, 0.85)
)

# Threshold lines in scaled coordinates (to match v2 guide lines)
x_cut_low <- scales::rescale(
  -fc_cut,
  from = range(volcano_data$log2fc),
  to = c(0.08, 0.92)
)
x_cut_high <- scales::rescale(
  fc_cut,
  from = range(volcano_data$log2fc),
  to = c(0.08, 0.92)
)
y_cut <- scales::rescale(
  p_cut,
  from = range(volcano_data$neg_log10_pval),
  to = c(0.10, 0.85)
)

volcano_data$color <- NA
volcano_data$color[volcano_data$regulation == "Not significant"] <- "#999999" # Gray
volcano_data$color[volcano_data$regulation == "Significant"] <- "#009E73" # Green
volcano_data$color[volcano_data$regulation == "Upregulated"] <- "#0072B2" # Blue
volcano_data$color[volcano_data$regulation == "Downregulated"] <- "#E69F00" # Orange/Yellow

# ---- Heatmap layer (v3) ----
rows <- paste0("G", sprintf("%02d", 1:12))
cols <- paste0("S", sprintf("%02d", 1:10))

base_matrix <- matrix(rnorm(12 * 10, mean = 0, sd = 0.5), nrow = 12, ncol = 10)
base_matrix[1:4, 1:5] <- base_matrix[1:4, 1:5] + 2.0 # up cluster
base_matrix[9:12, 6:10] <- base_matrix[9:12, 6:10] - 2.0 # down cluster
base_matrix[5:8, 4:7] <- base_matrix[5:8, 4:7] + 1.0 # moderate cluster

heatmap_df <- expand.grid(Row = rows, Col = cols)
heatmap_df$Value <- as.vector(base_matrix)
heatmap_df$Value <- scales::rescale(heatmap_df$Value, to = c(-2, 2))

heat_colors <- c("#0072B2", "#F7F7F7", "#E69F00")

# Map heatmap to the same 0-1 space as the volcano plot
heatmap_df$x <- scales::rescale(
  as.numeric(factor(heatmap_df$Col, levels = cols)),
  to = c(0.10, 0.90)
)
heatmap_df$y <- scales::rescale(
  as.numeric(factor(heatmap_df$Row, levels = rows)),
  to = c(0.16, 0.96)
)
tile_width <- (0.90 - 0.10) / length(cols)
tile_height <- (0.96 - 0.16) / length(rows)

# ---- Combined plot ----
# Heatmap as background, volcano points as foreground
p <- ggplot() +
  geom_tile(
    data = heatmap_df,
    aes(x = x, y = y, fill = Value),
    width = tile_width,
    height = tile_height,
    color = "#D0D0D0",
    linewidth = 0.2,
    alpha = 0.4,
    show.legend = FALSE
  ) +
  scale_fill_gradient2(
    low = heat_colors[1],
    mid = heat_colors[2],
    high = heat_colors[3],
    midpoint = 0,
    guide = "none"
  ) +
  # Subtle cow emoji watermark
  emoji_layer +
  # Volcano guide lines (from v2)
  geom_vline(
    xintercept = c(x_cut_low, x_cut_high),
    color = "#34495E",
    linewidth = 0.4,
    linetype = "dashed",
    alpha = 0.2
  ) +
  geom_hline(
    yintercept = y_cut,
    color = "#34495E",
    linewidth = 0.4,
    linetype = "dashed",
    alpha = 0.2
  ) +
  geom_vline(
    xintercept = 0.5,
    color = "#2C3E50",
    linewidth = 0.6,
    alpha = 0.2
  ) +
  # Volcano points on top
  geom_point(
    data = volcano_data,
    aes(x = x, y = y, color = color),
    size = 2.0,
    alpha = 0.17,
    shape = 16,
    show.legend = FALSE
  ) +
  scale_color_identity() +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = NA, color = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

# Create the hex sticker
sticker(
  p,
  package = "MOSuite",
  p_size = 24,
  p_color = "#1E3A8A",
  p_family = "sans",
  p_fontface = "bold",
  s_x = 1.0,
  s_y = 0.95,
  s_width = 2,
  s_height = 2,
  h_size = 1.3,
  h_color = "#1E90FF",
  h_fill = "#F8FBFF",
  dpi = 300,
  filename = "./inst/extdata/logo/logo_raw.png"
)
