# plot pairwise metrics heatmaps
# generates a lower-triangular grid of hex heatmaps colored by mean SND or SFD
# can be configured for 3 or 4 species

# setup ------------------------------------------------------------------

# load setup
source("code/01-setup.R")

# pick number of species (3 or 4)
richness <- 3

# pick target metric ("SND" or "SFD")
target_metric <- "SND"

# parameters
domi <- "IIB"
ratio <- "PNB"

# set metric-specific aesthetics
if (target_metric == "SND") {
  color_palette <- "A" # magma
  legend_title <- "Structural niche\ndifferences (SND)"
} else if (target_metric == "SFD") {
  color_palette <- "D" # viridis
  legend_title <- "Structural fitness\ndifferences (SFD)"
} else {
  stop("target_metric must be either 'SND' or 'SFD'")
}

# load empirical dataset
df <- read.table(
  paste0(
    "data/processed/empirical/str-coex-modules-",
    richness,
    "-species.txt"
  ),
  header = TRUE,
  sep = "\t"
)


# plotting function ------------------------------------------------------

# create a reusable function to generate hex heatmaps
plot_hex_heatmap <- function(data, x_col, y_col, x_lab, y_lab, metric, pal) {
  ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_vline(xintercept = 0, lty = "dotted", color = "grey30") +
    geom_hline(yintercept = 0, lty = "dotted", color = "grey30") +
    stat_summary_hex(
      # dynamically evaluate the target metric column
      aes(z = .data[[metric]]),
      fun = mean,
      bins = 40,
      color = "white",
      linewidth = 0.1
    ) +
    scale_fill_viridis_c(
      option = pal,
      direction = -1,
      begin = 0.2,
      end = 0.95,
      limits = c(0, round(max(data[[metric]], na.rm = TRUE)))
    ) +
    labs(x = x_lab, y = y_lab) +
    scale_x_continuous(position = "top") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 14, face = "bold", color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      panel.border = element_blank()
    )
}


# generate individual panels ---------------------------------------------

# row 1 (y = IIB)
p_pnb_iib <- plot_hex_heatmap(
  df,
  "pnb",
  "iib",
  ratio,
  domi,
  target_metric,
  color_palette
)
p_skew_iib <- plot_hex_heatmap(
  df,
  "skewness",
  "iib",
  "Skewness",
  domi,
  target_metric,
  color_palette
)
p_kurt_iib <- plot_hex_heatmap(
  df,
  "kurtosis",
  "iib",
  "Kurtosis",
  domi,
  target_metric,
  color_palette
)

# row 2 (y = Kurtosis)
p_pnb_kurt <- plot_hex_heatmap(
  df,
  "pnb",
  "kurtosis",
  ratio,
  "Kurtosis",
  target_metric,
  color_palette
)
p_skew_kurt <- plot_hex_heatmap(
  df,
  "skewness",
  "kurtosis",
  "Skewness",
  "Kurtosis",
  target_metric,
  color_palette
)

# row 3 (y = Skewness)
p_pnb_skew <- plot_hex_heatmap(
  df,
  "pnb",
  "skewness",
  ratio,
  "Skewness",
  target_metric,
  color_palette
)


# axis removal modifiers -------------------------------------------------

# define theme modifiers to strip only the titles from the inner panels
no_x <- theme(axis.title.x = element_blank())
no_y <- theme(axis.title.y = element_blank())

# apply modifiers to strip inner axes titles
p_skew_iib <- p_skew_iib + no_y
p_kurt_iib <- p_kurt_iib + no_y
p_pnb_kurt <- p_pnb_kurt + no_x
p_skew_kurt <- p_skew_kurt + no_x + no_y
p_pnb_skew <- p_pnb_skew + no_x


# legend panel -----------------------------------------------------------

# force the legend into a plot panel
snd_leg_raster <- ggplot(df, aes(x = skewness, y = kurtosis)) +
  # dynamically map fill to the target metric
  geom_point(
    aes(fill = .data[[target_metric]]),
    alpha = 0,
    shape = 21,
    color = NA
  ) +
  scale_fill_viridis_c(
    option = color_palette,
    direction = -1,
    begin = 0.2,
    end = 0.95,
    limits = c(0, round(max(df[[target_metric]], na.rm = TRUE))),
    name = legend_title,
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(8, "lines"),
      barheight = unit(1, "lines")
    )
  ) +
  theme_void() +
  theme(
    legend.position = c(0.5, 0.5),
    legend.direction = "horizontal",
    legend.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      vjust = 1.5
    ),
    legend.text = element_text(size = 12)
  )


# assemble final figure --------------------------------------------------

# triangular layout
arranged_snd_pairs <- ggpubr::ggarrange(
  p_pnb_iib,
  p_skew_iib,
  p_kurt_iib,
  p_pnb_kurt,
  p_skew_kurt,
  snd_leg_raster,
  p_pnb_skew,
  ggplot() + theme_void(),
  ggplot() + theme_void(),
  align = "hv",
  ncol = 3,
  nrow = 3,
  labels = c("a", "b", "c", "d", "e", "", "f", "", ""),
  label.x = 0.055,
  label.y = 0.95,
  font.label = list(size = 18, face = "bold")
)


# export -----------------------------------------------------------------

# save the figure dynamically named by metric and species count
ggsave(
  filename = paste0(
    "suppl-mat/suppl-figures/fig-heatmaps-",
    target_metric,
    "-",
    richness,
    "-species.jpeg"
  ),
  plot = arranged_snd_pairs,
  device = "jpeg",
  dpi = dpi,
  height = 10,
  width = 10,
  bg = "white"
)

# info
cat(
  "Done! Generated hex heatmaps for",
  target_metric,
  "with",
  richness,
  "species.\n"
)
