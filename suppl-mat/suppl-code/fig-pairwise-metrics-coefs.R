# evaluate relationships between pairwise network metrics and structural outputs

# setup ------------------------------------------------------------------

# load setup
source("code/01-setup.R")

# define richness level (3 and 4)
richness_levels <- c(3, 4)

# define labels
metric_labels <- c(
  "iib" = "IIB",
  "pnb" = "PNB",
  "skewness" = "Skewness",
  "kurtosis" = "Kurtosis"
)

# define pairwise combinations
pair_list <- list(
  list(vars = c("iib", "pnb"), name = "IIB & PNB"),
  list(vars = c("iib", "skewness"), name = "IIB & Skewness"),
  list(vars = c("iib", "kurtosis"), name = "IIB & Kurtosis"),
  list(vars = c("pnb", "skewness"), name = "PNB & Skewness"),
  list(vars = c("pnb", "kurtosis"), name = "PNB & Kurtosis"),
  list(vars = c("skewness", "kurtosis"), name = "Skewness & Kurtosis")
)


# processing function ----------------------------------------------------

# function to process each richness level and extract coefs for pairwise models
process_pairwise <- function(richness) {
  str_coex <- read.table(
    paste0(
      "data/processed/empirical/str-coex-modules-",
      richness,
      "-species.txt"
    ),
    header = TRUE,
    sep = "\t"
  )

  # standardize predictors
  str_coex_scaled <- str_coex %>%
    mutate(
      iib = as.numeric(scale(iib)),
      pnb = as.numeric(scale(pnb)),
      skewness = as.numeric(scale(skewness)),
      kurtosis = as.numeric(scale(kurtosis))
    )

  # iterate over pairs and extract coefficients
  res_list <- lapply(pair_list, function(p) {
    f_snd <- as.formula(paste("SND ~", p$vars[1], "+", p$vars[2]))
    f_sfd <- as.formula(paste("SFD ~", p$vars[1], "+", p$vars[2]))

    m_snd <- lm(f_snd, data = str_coex_scaled)
    m_sfd <- lm(f_sfd, data = str_coex_scaled)

    d_snd <- tidy(m_snd, conf.int = TRUE, conf.level = 0.99) %>%
      mutate(response = "SND")
    d_sfd <- tidy(m_sfd, conf.int = TRUE, conf.level = 0.99) %>%
      mutate(response = "SFD")

    bind_rows(d_snd, d_sfd) %>%
      filter(term != "(Intercept)") %>%
      mutate(model_pair = p$name)
  })

  bind_rows(res_list) %>%
    mutate(module_size = paste(richness, "species")) %>%
    dplyr::select(
      response,
      model_pair,
      term,
      estimate,
      ci_low = conf.low,
      ci_high = conf.high,
      module_size
    )
}


# compute and extract data -----------------------------------------------

# bind results for all richness levels
all_results <- bind_rows(lapply(richness_levels, process_pairwise))

# prep factors for the plot
plot_data <- all_results %>%
  mutate(
    term = factor(term, levels = rev(c("iib", "pnb", "skewness", "kurtosis"))),
    module_size = factor(module_size, levels = c("4 species", "3 species")),
    model_pair = factor(model_pair, levels = sapply(pair_list, `[[`, "name"))
  )


# plotting ---------------------------------------------------------------

# shared base theme with updated aesthetics
base_theme <- theme_classic() +
  theme(
    # text to black
    text = element_text(color = "black"),

    # draw only bottom and left axes
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),

    # strip settings (removes facet titles and backgrounds)
    strip.background = element_blank(),
    strip.text = element_blank(),

    # text formatting
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title.x = element_text(color = "black", size = 12),

    # add vertical space between panels so the inner axes don't overlap
    panel.spacing = unit(1.2, "lines")
  )

# panel a: snd plot
plot_snd <- ggplot(
  plot_data %>% filter(response == "SND"),
  aes(x = estimate, y = term, fill = module_size)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(
    aes(xmin = ci_low, xmax = ci_high),
    shape = 21,
    color = "#be3872",
    size = 0.55,
    stroke = 0.8,
    position = position_dodge(width = 0.5)
  ) +
  scale_y_discrete(labels = metric_labels) +
  scale_fill_manual(
    values = c("3 species" = "white", "4 species" = "#be3872")
  ) +
  facet_wrap(
    ~model_pair,
    ncol = 1,
    scales = "free_y",
    axes = "all",
    axis.labels = "all"
  ) +
  guides(
    fill = guide_legend(
      title = "Module size",
      reverse = TRUE,
      override.aes = list(
        shape = 21,
        color = "black",
        fill = c("white", "black"),
        stroke = 0.8
      )
    )
  ) +
  labs(x = "Effect on structural niche differences\n(SND [99% CI])", y = NULL) +
  base_theme

# panel b: sfd plot
plot_sfd <- ggplot(
  plot_data %>% filter(response == "SFD"),
  aes(x = estimate, y = term, fill = module_size)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(
    aes(xmin = ci_low, xmax = ci_high),
    shape = 21,
    color = "#208f8a",
    size = 0.55,
    stroke = 0.8,
    position = position_dodge(width = 0.5)
  ) +
  scale_y_discrete(labels = metric_labels) +
  scale_fill_manual(
    values = c("3 species" = "white", "4 species" = "#208f8a")
  ) +
  facet_wrap(
    ~model_pair,
    ncol = 1,
    scales = "free_y",
    axes = "all",
    axis.labels = "all"
  ) +
  guides(fill = "none") +
  labs(
    x = "Effect on structural fitness differences\n(SFD [99% CI])",
    y = NULL
  ) +
  base_theme +
  theme(
    axis.text.y = element_blank()
  )

# combine using patchwork
plot_combined <- plot_snd +
  plot_sfd +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 16, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12, color = "black"),
    legend.text = element_text(size = 12, color = "black")
  )


# export -----------------------------------------------------------------

# save the combined plot as a high-resolution JPEG
ggsave(
  filename = "suppl-mat/suppl-figures/fig-pairwise-metrics-coefs.jpeg",
  plot = plot_combined,
  device = "jpeg",
  dpi = dpi,
  height = 9.5,
  width = 8,
  bg = "white"
)

# info
cat("Done! Pairwise coefficient forest plot saved.\n")
