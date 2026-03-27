# evaluate relationships between network metrics and structural outputs

# setup ------------------------------------------------------------------

# load setup and libraries
source("code/01-setup.R")

# define module sizes to process
richness_levels <- c(3, 4)

# define metric labels for plotting
metric_labels <- c(
  "iid" = "IID",
  "pnd" = "PND",
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
      iid = as.numeric(scale(iid)),
      pnd = as.numeric(scale(pnd)),
      skewness = as.numeric(scale(skewness)),
      kurtosis = as.numeric(scale(kurtosis))
    )

  # run glms
  m_snd <- glm(SND ~ iid + pnd + skewness + kurtosis, data = str_coex_scaled)
  m_sfd <- glm(SFD ~ iid + pnd + skewness + kurtosis, data = str_coex_scaled)

  # extract cis and tag with richness
  bind_rows(
    tidy(m_snd, conf.int = TRUE, conf.level = 0.99) %>%
      mutate(response = "SND"),
    tidy(m_sfd, conf.int = TRUE, conf.level = 0.99) %>% mutate(response = "SFD")
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
    term = factor(term, levels = rev(c("iid", "pnd", "skewness", "kurtosis"))),
    # Ensure proper ordering so 3 species dodges "above" 4 species visually
    module_size = factor(module_size, levels = c("4 species", "3 species"))
  )


# plotting ---------------------------------------------------------------

# define the base theme explicitly so it never fails
base_theme <- theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

# panel a: snd plot
plot_snd <- ggplot(
  plot_data %>% filter(response == "SND"),
  aes(x = estimate, y = term, fill = module_size, group = module_size)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(
    aes(xmin = ci_low, xmax = ci_high),
    shape = 21,
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
  guides(
    fill = guide_legend(
      title = "Module size",
      override.aes = list(
        shape = 21,
        color = "black",
        fill = c("white", "black"),
        stroke = 0.8
      )
    )
  ) +
  labs(
    x = "Effect on structural niche differences\n(SND [99% CI])",
    y = NULL
  ) +
  # add the letter a inside the top left of the panel
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
  base_theme

# panel b: sfd plot
plot_sfd <- ggplot(
  plot_data %>% filter(response == "SFD"),
  aes(x = estimate, y = term, fill = module_size, group = module_size)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(
    aes(xmin = ci_low, xmax = ci_high),
    shape = 21,
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
  guides(fill = "none") +
  labs(
    x = "Effect on structural fitness differences\n(SFD [99% CI])",
    y = NULL
  ) +
  # add the letter b inside the top left of the panel
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
    axis.text.y = element_blank()
  )

# combine using patchwork
plot_combined <- plot_snd +
  plot_sfd +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )


# export -----------------------------------------------------------------

ggsave(
  filename = "results/figures/fig-glm-all-additive-effects.jpeg",
  plot = plot_combined,
  device = "jpeg",
  dpi = dpi,
  height = 4,
  width = 6.5
)

# info
cat(
  "Done!\n"
)
