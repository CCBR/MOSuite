mosuite_palette <- c(
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

usethis::use_data(mosuite_palette, overwrite = TRUE)

ggplot2::ggsave(
  filename = "mosuite_palette.png",
  plot = display_palette(mosuite_palette),
  path = "man/figures",
  width = 7,
  height = 2,
  units = "in",
  dpi = 300
)
