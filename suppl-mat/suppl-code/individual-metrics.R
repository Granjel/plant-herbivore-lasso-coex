# evaluate individual relationships between metrics and structural outputs

# setup ------------------------------------------------------------------

# load setup
source("code/01-setup.R")

# define richness levels and metrics to analyse
richness_levels <- c(3, 4)
metrics <- c("iid", "pnd", "skewness", "kurtosis")
metric_labels <- c(
  "iid" = "IID",
  "pnd" = "PND",
  "skewness" = "Skewness",
  "kurtosis" = "Kurtosis"
)


# processing function ----------------------------------------------------

# fn to fit a single model for each metric and response, extract, and combine
process_individual_models <- function(richness) {
  str_coex <- read.table(
    paste0(
      "data/processed/empirical/str-coex-modules-",
      richness,
      "-species.txt"
    ),
    header = TRUE,
    sep = "\t"
  )

  str_coex_scaled <- str_coex %>%
    mutate(across(all_of(metrics), ~ as.numeric(scale(.x))))

  lapply(metrics, function(m) {
    m_snd <- glm(as.formula(paste("SND ~", m)), data = str_coex_scaled)
    m_sfd <- glm(as.formula(paste("SFD ~", m)), data = str_coex_scaled)

    bind_rows(
      tidy(m_snd, conf.int = TRUE, conf.level = 0.99) %>%
        mutate(response = "SND"),
      tidy(m_sfd, conf.int = TRUE, conf.level = 0.99) %>%
        mutate(response = "SFD")
    ) %>%
      mutate(metric = m, module_size = paste(richness, "species"))
  }) %>%
    bind_rows()
}


# compute data -----------------------------------------------------------

all_results <- bind_rows(lapply(richness_levels, process_individual_models))

# isolate ONLY intercepts for the supplementary table
intercept_table <- all_results %>%
  filter(term == "(Intercept)") %>%
  dplyr::select(module_size, metric, response, estimate, conf.low, conf.high)

write.table(
  intercept_table,
  file = "suppl-mat/suppl-tables/table-s-individual-glm-intercepts.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# prep plot data
plot_data <- all_results %>%
  filter(term != "(Intercept)") %>%
  mutate(
    metric = factor(metric, levels = metrics),
    module_size = factor(module_size, levels = c("4 species", "3 species"))
  )


# plotting ---------------------------------------------------------------

# helper function for the individual "strips"
make_mini_plot <- function(
  m_name,
  resp_name,
  color_hex,
  show_y = FALSE,
  show_x = FALSE,
  tag = NULL
) {
  df_sub <- plot_data %>% filter(metric == m_name, response == resp_name)

  x_label <- if (show_x) {
    if (resp_name == "SND") {
      "Effect on structural niche differences\n(SND [99% CI])"
    } else {
      "Effect on structural fitness differences\n(SFD [99% CI])"
    }
  } else {
    NULL
  }

  ggplot(df_sub, aes(x = estimate, y = metric, fill = module_size)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_pointrange(
      aes(xmin = conf.low, xmax = conf.high),
      shape = 21,
      color = color_hex,
      size = 0.55,
      stroke = 0.8,
      position = position_dodge(width = 0.4)
    ) +
    scale_y_discrete(labels = metric_labels) +
    scale_fill_manual(
      values = c("3 species" = "white", "4 species" = color_hex),
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
    labs(x = x_label, y = NULL, title = tag) +
    theme_classic() +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = -0.15,
        vjust = 1
      ),
      axis.text.y = if (show_y) {
        element_text(size = 10, color = "black")
      } else {
        element_blank()
      },
      axis.line.y = if (show_y) element_line() else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank(),
      axis.text.x = element_text(size = 8),
      axis.title.x = if (show_x) {
        element_text(size = 10, face = "plain")
      } else {
        element_blank()
      },
      plot.margin = margin(10, 5, 10, 5)
    )
}

# row tags
row_tags <- c("a", "b", "c", "d")

# generate the 8 plots - tags applied only to the SND (left) column
plots_snd <- lapply(1:4, function(i) {
  make_mini_plot(
    metrics[i],
    "SND",
    "#be3872",
    show_y = TRUE,
    show_x = (i == 4),
    tag = row_tags[i]
  )
})

plots_sfd <- lapply(1:4, function(i) {
  make_mini_plot(
    metrics[i],
    "SFD",
    "#208f8a",
    show_y = FALSE,
    show_x = (i == 4),
    tag = NULL
  )
})

# assemble row by row
row1 <- plots_snd[[1]] | plots_sfd[[1]]
row2 <- plots_snd[[2]] | plots_sfd[[2]]
row3 <- plots_snd[[3]] | plots_sfd[[3]]
row4 <- plots_snd[[4]] | plots_sfd[[4]]

# combine all rows into the final grid
plot_final <- (row1 / row2 / row3 / row4) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )


# export -----------------------------------------------------------------

# save the plot
ggsave(
  filename = "suppl-mat/suppl-figures/fig-individual-metrics.jpeg",
  plot = plot_final,
  device = "jpeg",
  dpi = 320,
  height = 7,
  width = 6.25
)

# info
cat("done! row-tagged supplementary grid generated.\n")
