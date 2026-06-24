# evaluate relationships between network metrics and structural outputs

# setup ------------------------------------------------------------------

# load setup and libraries
source("code/01-setup.R")
library(cowplot)

# define module sizes to process
richness_levels <- c(3, 4)

# define metric labels for plotting
metric_labels <- c(
  "iib" = "Intra-Interspecific\nBalance (IIB)",
  "pnb" = "Positive-Negative\nBalance (PNB)",
  "skewness" = "Skewness",
  "kurtosis" = "Kurtosis"
)

# processing function ----------------------------------------------------

# function to load data, run models, and extract theoretical 99% cis
process_modules <- function(richness) {
  # load empirical structural results
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

  # run glms using tweedie distribution
  m_snd <- glmmTMB(
    SND ~ iib + pnb + skewness + kurtosis,
    family = tweedie(link = "log"),
    data = str_coex_scaled
  )

  m_sfd <- glmmTMB(
    SFD ~ iib + pnb + skewness + kurtosis,
    family = tweedie(link = "log"),
    data = str_coex_scaled
  )

  # check model diagnostics
  # simulateResiduals(fittedModel = m_snd, n = 250, plot = TRUE)
  # simulateResiduals(fittedModel = m_sfd, n = 250, plot = TRUE)

  # extract cis and tag with richness
  bind_rows(
    tidy(m_snd, conf.int = TRUE, conf.level = 0.99, effects = "fixed") %>%
      mutate(response = "SND"),
    tidy(m_sfd, conf.int = TRUE, conf.level = 0.99, effects = "fixed") %>%
      mutate(response = "SFD")
  ) %>%
    dplyr::select(
      response,
      term,
      estimate,
      ci_low = conf.low,
      ci_high = conf.high
    ) %>%
    mutate(module_size = paste(richness, "species"))
}

# compute and extract data -----------------------------------------------

# map function across both richness levels and bind together
all_results <- bind_rows(lapply(richness_levels, process_modules))

# isolate intercepts and save to a table
intercept_table <- all_results %>%
  filter(term == "(Intercept)") %>%
  dplyr::select(module_size, response, estimate, ci_low, ci_high)

# export the intercepts to a table for reporting
write.table(
  intercept_table,
  file = "results/tables/table-glm-all-additive-intercepts.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# filter out intercepts and prep factors for the plot
plot_data <- all_results %>%
  filter(term != "(Intercept)") %>%
  mutate(
    hypothesis_support = case_when(
      # expected effects for snd are positive
      response == "SND" & estimate > 0 ~ "Yes",
      response == "SND" & estimate <= 0 ~ "No",

      # expected effects for sfd are negative for iib and pnb
      response == "SFD" & term %in% c("iib", "pnb") & estimate < 0 ~ "Yes",
      response == "SFD" & term %in% c("iib", "pnb") & estimate >= 0 ~ "No",

      # unclear expectation for skewness and kurtosis effects on sfd
      response == "SFD" & term %in% c("skewness", "kurtosis") ~
        "Unclear expectation"
    ),
    hypothesis_support = factor(
      hypothesis_support,
      levels = c("Yes", "No", "Unclear expectation")
    ),
    term = factor(
      term,
      levels = rev(c("iib", "pnb", "skewness", "kurtosis"))
    ),
    # ensure proper ordering so 3 species dodges above 4 species visually
    module_size = factor(module_size, levels = c("4 species", "3 species"))
  )

# plotting ---------------------------------------------------------------

# define the base theme explicitly so it never fails
base_theme <- theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    legend.box = "horizontal"
  )

# panel a: snd plot
plot_snd <- ggplot(
  plot_data %>% filter(response == "SND"),
  aes(
    x = estimate,
    y = term,
    fill = module_size,
    shape = hypothesis_support,
    group = module_size
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_pointrange(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    color = "#be3872",
    size = 0.55,
    stroke = 0.8,
    position = position_dodge(width = 0.4)
  ) +
  scale_y_discrete(labels = metric_labels) +
  scale_fill_manual(
    values = c("3 species" = "white", "4 species" = "#be3872"),
    breaks = c("3 species", "4 species")
  ) +
  scale_shape_manual(
    values = c(
      "Yes" = 21,
      "No" = 22,
      "Unclear expectation" = 24
    ),
    breaks = c("Yes", "No", "Unclear expectation"),
    drop = FALSE
  ) +
  guides(
    fill = guide_legend(
      title = "Module size",
      override.aes = list(
        shape = 21,
        color = "black",
        fill = c("white", "black"),
        stroke = 0.8
      )
    ),
    shape = "none"
  ) +
  labs(
    x = "Association with structural\nniche differences",
    y = NULL
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "a",
    fontface = "bold",
    size = 5,
    hjust = -1,
    vjust = 1.25
  ) +
  base_theme +
  theme(
    legend.position = "bottom"
  )

# panel b: sfd plot
plot_sfd <- ggplot(
  plot_data %>% filter(response == "SFD"),
  aes(
    x = estimate,
    y = term,
    fill = module_size,
    shape = hypothesis_support,
    group = module_size
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_pointrange(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    color = "#208f8a",
    size = 0.55,
    stroke = 0.8,
    position = position_dodge(width = 0.4)
  ) +
  scale_y_discrete(labels = metric_labels) +
  scale_fill_manual(
    values = c("3 species" = "white", "4 species" = "#208f8a"),
    breaks = c("3 species", "4 species")
  ) +
  scale_shape_manual(
    values = c(
      "Yes" = 21,
      "No" = 22,
      "Unclear expectation" = 24
    ),
    breaks = c("Yes", "No", "Unclear expectation"),
    drop = FALSE
  ) +
  guides(
    fill = "none",
    shape = guide_legend(
      title = "Matches Fig. 2 expectations?",
      override.aes = list(
        fill = "black",
        color = "black",
        stroke = 0.8,
        linetype = 0,
        linewidth = 0
      )
    )
  ) +
  labs(
    x = "Association with structural\nfitness differences",
    y = NULL
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "b",
    fontface = "bold",
    size = 5,
    hjust = -1,
    vjust = 1.25
  ) +
  base_theme +
  theme(
    axis.text.y = element_blank(),
    legend.position = "top"
  )

# extract legends so they can span the whole figure
legend_top <- cowplot::get_legend(plot_sfd)
legend_bottom <- cowplot::get_legend(plot_snd)

# remove legends from the actual panels
plot_snd_panel <- plot_snd + theme(legend.position = "none")
plot_sfd_panel <- plot_sfd + theme(legend.position = "none")

# combine panels
plots_row <- plot_snd_panel +
  plot_sfd_panel +
  plot_layout(widths = c(1, 1))

# assemble full figure with full-width legends
plot_combined <- (patchwork::wrap_elements(legend_top) /
  plots_row /
  patchwork::wrap_elements(legend_bottom)) +
  plot_layout(heights = c(0.16, 1, 0.14))

# export -----------------------------------------------------------------

ggsave(
  filename = "results/figures/fig-glm-all-additive.jpeg",
  plot = plot_combined,
  device = "jpeg",
  dpi = dpi,
  height = 4.3,
  width = 6.8
)

# info
cat("Done!\n")
