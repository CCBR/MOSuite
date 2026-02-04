#!/usr/bin/env Rscript
# Create hex logo for MOSuite package
# Design: Combine volcano plot + heatmap

library(ggplot2)
library(hexSticker)
library(ggimage)
library(magick)
library(tidyverse)

# Create directory if it doesn't exist

set.seed(123)

emoji_url <- "https://raw.githubusercontent.com/twitter/twemoji/master/assets/72x72/1f42e.png"
emoji_temp <- tempfile(fileext = ".png")
download.file(emoji_url, emoji_temp, quiet = TRUE, mode = "wb")

# Convert to grayscale and apply transparency
emoji_gray <- magick::image_read(emoji_temp) |>
  magick::image_quantize(colorspace = "gray")
emoji_gray_transparent <- magick::image_fx(
  emoji_gray,
  "0.3*u",
  channel = "alpha"
)
magick::image_write(emoji_gray_transparent, emoji_temp)

emoji_df <- data.frame(
  x = 0.366,
  y = 0.73,
  image = emoji_temp
)


# ---- Volcano plot layer ----

load(here::here('data', 'nidap_deg_analysis.rda'))
volcano_data <- nidap_deg_analysis |>
  rename(log2fc = `C-A_logFC`) |>
  mutate(neg_log10_pval = -log10(`C-A_pval`))

fc_cut <- 1.0
p_cut <- 2.5

volcano_data <- volcano_data |>
  mutate(
    regulation = case_when(
      log2fc >= fc_cut & neg_log10_pval >= p_cut ~ "Upregulated",
      log2fc <= -fc_cut & neg_log10_pval >= p_cut ~ "Downregulated",
      neg_log10_pval >= p_cut ~ "Significant",
      .default = "Not significant"
    ),
    x = scales::rescale(log2fc, to = c(0.12, 0.88)),
    y = scales::rescale(neg_log10_pval, to = c(0.18, 0.94)),
    color = case_when(
      regulation == "Upregulated" ~ "#4e9db5",
      regulation == "Downregulated" ~ "#ecba4c",
      regulation == "Significant" ~ "#528230",
      regulation == "Not significant" ~ "#999999",
      .default = NA_character_
    )
  )

# ---- Heatmap layer ----
rows <- paste0("G", sprintf("%02d", 1:12))
cols <- paste0("S", sprintf("%02d", 1:10))

base_matrix <- matrix(rnorm(12 * 10, mean = 0, sd = 0.5), nrow = 12, ncol = 10)
base_matrix[1:4, 1:5] <- base_matrix[1:4, 1:5] + 2.0 # up cluster
base_matrix[9:12, 6:10] <- base_matrix[9:12, 6:10] - 2.0 # down cluster
base_matrix[5:8, 4:7] <- base_matrix[5:8, 4:7] + 1.0 # moderate cluster

heatmap_df <- expand.grid(Row = rows, Col = cols) |>
  mutate(
    Value = as.vector(base_matrix),
    Value = scales::rescale(Value, to = c(-2, 2)),
    x = scales::rescale(
      as.numeric(factor(Col, levels = cols)),
      to = c(0.10, 0.90)
    ),
    y = scales::rescale(
      as.numeric(factor(Row, levels = rows)),
      to = c(0.16, 0.96)
    )
  )

heat_colors <- c("#4e9db5", "#F7F7F7", "#ecba4c")
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
    alpha = 0.5,
    show.legend = FALSE
  ) +
  scale_fill_gradient2(
    low = heat_colors[1],
    mid = heat_colors[2],
    high = heat_colors[3],
    midpoint = 0,
    guide = "none"
  ) +
  # Volcano points on top
  geom_point(
    data = volcano_data,
    aes(x = x, y = y, color = color),
    size = 2.0,
    alpha = 0.5,
    shape = 16,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = 0.5,
    y = 0.725,
    label = "MOSuite",
    size = 16,
    fontface = "bold",
    family = "sans",
    color = "#296b7f"
  ) +
  ggimage::geom_image(
    # emoji
    data = emoji_df,
    aes(x = 0.33, y = 0.725, image = image),
    size = 0.08,
    inherit.aes = FALSE
  ) +
  scale_color_identity() +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "#FFF", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  )
print(p)

ggsave(
  here::here('inst', 'extdata', 'logo', 'mosuite_logo_with_text.png'),
  plot = p,
  dpi = 300,
  width = 3,
  height = 3
)

# background only
p_background <- ggplot() +
  geom_tile(
    data = heatmap_df,
    aes(x = x, y = y, fill = Value),
    width = tile_width,
    height = tile_height,
    color = "#D0D0D0",
    linewidth = 0.2,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  scale_fill_gradient2(
    low = heat_colors[1],
    mid = heat_colors[2],
    high = heat_colors[3],
    midpoint = 0,
    guide = "none"
  ) +
  # Volcano points on top
  geom_point(
    data = volcano_data,
    aes(x = x, y = y, color = color),
    size = 2.0,
    alpha = 0.7,
    shape = 16,
    show.legend = FALSE
  ) +
  ggimage::geom_image(
    # emoji
    data = emoji_df,
    aes(x = 0.33, y = 0.725, image = image),
    size = 0.08,
    inherit.aes = FALSE
  ) +
  scale_color_identity() +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "#FFF", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  )
print(p_background)
ggsave(
  here::here('inst', 'extdata', 'logo', 'mosuite_logo_background.png'),
  plot = p_background,
  dpi = 300,
  width = 3,
  height = 3
)

# need to run this line interactively for it to actually overwrite the logo file
usethis::use_logo(here::here(
  'inst',
  'extdata',
  'logo',
  'mosuite_logo_background.png'
))
pkgdown::build_favicons(overwrite = TRUE)
