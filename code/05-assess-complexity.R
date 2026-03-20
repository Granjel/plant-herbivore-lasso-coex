# plot empirical and bootstrapped interactions (figure 3)

# setup ------------------------------------------------------------------

library(tidyverse)
library(ggpubr)

# load empirical and bootstrapped matrices
load("data/processed/empirical/empirical-matrices.RData")
load("data/processed/bootstrapped/bootstrapped-matrices.RData")

# define ci to crop plot and regression lines color
perc <- 0.01
color_smooth <- c("darkgreen", "deepskyblue4")[1]
transparency <- 0.15

# abbreviations (ensure this matches your 'species' vector order)
species_abbrev <- c(
  "Am",
  "Ao",
  "Ae",
  "Be",
  "Cj",
  "Ca",
  "Cr",
  "Dg",
  "Dc",
  "Er",
  "Em",
  "Fa",
  "Fr",
  "Gv",
  "Gd",
  "Gr",
  "Lv",
  "Lp",
  "Lc",
  "Ma",
  "Or",
  "Pe",
  "Ph",
  "Pl",
  "Pa",
  "Pp",
  "Pt",
  "Ra",
  "Rx",
  "Sp",
  "So",
  "To",
  "Tf",
  "Tp",
  "Vo",
  "Vp"
)


# data preparation: igr and gamma ----------------------------------------

# 1. extract empirical igr and gamma
emp_igr_gamma <- tibble(
  spp_code = species,
  spp = species_abbrev,
  bootstrapped = FALSE,
  igr = igr,
  grasshoppers = rowSums(gamma)
)

# 2. extract bootstrapped igr and gamma
boot_igr_gamma <- boot_matrices %>%
  select(boot_id, igr, gamma) %>%
  mutate(
    spp_code = list(species),
    grasshoppers = map(gamma, rowSums)
  ) %>%
  select(-gamma) %>%
  unnest(c(spp_code, igr, grasshoppers)) %>%
  mutate(bootstrapped = TRUE)

# 3. combine for scatter plot (panel a)
igr_gamma_for_plotting <- bind_rows(
  emp_igr_gamma %>% select(-spp),
  boot_igr_gamma %>% select(-boot_id)
)

# 4. calculate 95% CIs for forest plot (panel b)
igr_cis <- boot_igr_gamma %>%
  group_by(spp_code) %>%
  summarise(
    igr_low = quantile(igr, 0.025),
    igr_high = quantile(igr, 0.975),
    igr_gamma_low = quantile(igr + grasshoppers, 0.025),
    igr_gamma_high = quantile(igr + grasshoppers, 0.975),
    .groups = "drop"
  )

# 5. format data for forest plot (panel b) exactly as your original ggplot expects
igr_forest_data <- emp_igr_gamma %>%
  left_join(igr_cis, by = "spp_code") %>%
  # duplicate rows: one for plants alone, one for plants + gamma
  uncount(2, .id = "condition") %>%
  mutate(
    added_gamma = ifelse(condition == 1, FALSE, TRUE),
    plot_igr = ifelse(added_gamma, igr + grasshoppers, igr),
    plot_low = ifelse(added_gamma, igr_gamma_low, igr_low),
    plot_high = ifelse(added_gamma, igr_gamma_high, igr_high)
  )

# re-scale variables
resc <- max(igr_forest_data$plot_high)
igr_gamma_for_plotting$igr <- igr_gamma_for_plotting$igr / resc
igr_gamma_for_plotting$grasshoppers <- igr_gamma_for_plotting$grasshoppers /
  resc


# data preparation: alpha and beta ---------------------------------------

# 1. extract empirical alpha and beta
emp_alpha_beta <- tibble(
  bootstrapped = FALSE,
  alpha_coefficients = as.vector(alpha),
  beta_coefficients = as.vector(Reduce("+", beta_gp)) +
    as.vector(Reduce("+", beta_pp))
)

# 2. extract bootstrapped alpha and beta
boot_alpha_beta <- boot_matrices %>%
  mutate(
    alpha_coefficients = map(alpha, as.vector),
    beta_coefficients = map2(
      beta_pp,
      beta_gp,
      ~ as.vector(Reduce("+", .x)) + as.vector(Reduce("+", .y))
    )
  ) %>%
  select(alpha_coefficients, beta_coefficients) %>%
  unnest(c(alpha_coefficients, beta_coefficients)) %>%
  mutate(bootstrapped = TRUE)

# 3. combine for panels c and d
alpha_beta_for_plotting <- bind_rows(emp_alpha_beta, boot_alpha_beta) %>%
  mutate(changed = alpha_coefficients + beta_coefficients)


# plotting ---------------------------------------------------------------

# panel a: igr vs gamma
plot_igr_gamma <- ggplot(
  data = igr_gamma_for_plotting %>% filter(bootstrapped == TRUE),
  aes(x = igr, y = grasshoppers)
) +
  geom_point(alpha = transparency, color = "grey") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_x_continuous(
    limits = c(
      quantile(
        (igr_gamma_for_plotting %>% filter(bootstrapped == TRUE))$igr,
        perc / 2
      ),
      round(quantile(
        (igr_gamma_for_plotting %>% filter(bootstrapped == TRUE))$igr,
        1 - (perc / 2)
      ))
    )
  ) +
  scale_y_continuous(
    limits = c(
      round(
        quantile(
          (igr_gamma_for_plotting %>%
            filter(bootstrapped == TRUE))$grasshoppers,
          perc / 2
        ),
        3
      ),
      round(
        quantile(
          (igr_gamma_for_plotting %>%
            filter(bootstrapped == TRUE))$grasshoppers,
          1 - (perc / 2)
        ),
        3
      )
    ),
    breaks = seq(
      round(
        quantile(
          (igr_gamma_for_plotting %>%
            filter(bootstrapped == TRUE))$grasshoppers,
          perc / 2
        ),
        2
      ),
      round(
        quantile(
          (igr_gamma_for_plotting %>%
            filter(bootstrapped == TRUE))$grasshoppers,
          1 - (perc / 2)
        ),
        2
      ),
      by = 0.02
    )
  ) +
  geom_point(
    data = igr_gamma_for_plotting %>% filter(bootstrapped == FALSE),
    aes(x = igr, y = grasshoppers),
    color = "black",
    shape = 19
  ) +
  geom_smooth(
    data = igr_gamma_for_plotting %>% filter(bootstrapped == TRUE),
    aes(x = igr, y = grasshoppers),
    method = "glm",
    color = color_smooth,
    linetype = "dashed",
    se = FALSE
  ) +
  geom_smooth(
    data = igr_gamma_for_plotting %>% filter(bootstrapped == FALSE),
    aes(x = igr, y = grasshoppers),
    method = "glm",
    color = color_smooth,
    se = FALSE
  ) +
  xlab("Intrinsic growth rate") +
  ylab(expression("Grasshoppers' effects (" * gamma * ")")) +
  theme_classic() +
  theme(axis.title.x = element_text(colour = "white"))


# panel b: igr forest plot
plot_igr_gamma_changed <- ggplot(
  igr_forest_data,
  aes(
    y = reorder(spp, -plot_igr, FUN = median),
    x = plot_igr / resc,
    color = added_gamma
  )
) +
  geom_pointrange(
    aes(xmin = plot_low / resc, xmax = plot_high / resc),
    position = position_dodge(width = -0.75),
    size = 0.04
  ) +
  scale_color_manual(
    name = "Value (95% CI)",
    values = c("black", "grey70"),
    labels = c(
      "Plants alone",
      expression("Plants + grasshoppers (" * gamma * ")")
    )
  ) +
  xlab("Intrinsic growth rate") +
  ylab("Plant species") +
  theme_classic() +
  theme(
    axis.text.x = element_text(hjust = 0.5, vjust = 1, angle = 0),
    legend.position = c(0.65, 0.875),
    axis.text.y = element_text(size = 7)
  )


# arrange panels a and b
arranged_igr <- ggpubr::ggarrange(
  plot_igr_gamma,
  plot_igr_gamma_changed,
  heights = c(0.4, 0.6),
  ncol = 1,
  nrow = 2,
  align = "v",
  labels = c("a", "b"),
  vjust = 0.8
)


# panel c: alpha vs beta
plot_alpha_beta <- ggplot(
  data = alpha_beta_for_plotting %>% filter(bootstrapped == TRUE),
  aes(x = alpha_coefficients, y = beta_coefficients)
) +
  geom_point(alpha = transparency, color = "grey") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(
    limits = c(
      quantile(
        (alpha_beta_for_plotting %>%
          filter(bootstrapped == TRUE))$alpha_coefficients,
        perc / 2
      ),
      quantile(
        (alpha_beta_for_plotting %>%
          filter(bootstrapped == TRUE))$alpha_coefficients,
        1 - (perc / 2)
      )
    )
  ) +
  scale_y_continuous(
    limits = c(
      quantile(
        (alpha_beta_for_plotting %>%
          filter(bootstrapped == TRUE))$beta_coefficients,
        perc / 2
      ),
      quantile(
        (alpha_beta_for_plotting %>%
          filter(bootstrapped == TRUE))$beta_coefficients,
        1 - (perc / 2)
      )
    )
  ) +
  geom_point(
    data = alpha_beta_for_plotting %>% filter(bootstrapped == FALSE),
    aes(x = alpha_coefficients, y = beta_coefficients),
    color = "black",
    shape = 19
  ) +
  geom_smooth(
    data = alpha_beta_for_plotting %>% filter(bootstrapped == TRUE),
    aes(x = alpha_coefficients, y = beta_coefficients),
    method = "glm",
    color = color_smooth,
    linetype = "dashed",
    se = FALSE
  ) +
  geom_smooth(
    data = alpha_beta_for_plotting %>% filter(bootstrapped == FALSE),
    aes(x = alpha_coefficients, y = beta_coefficients),
    method = "glm",
    color = color_smooth,
    se = FALSE
  ) +
  xlab(expression("Plant-plant interactions (" * alpha * ")")) +
  ylab(expression("Higher-order interactions (" * beta * ")")) +
  theme_classic() +
  theme(axis.title.x = element_text(color = "white"))


# panel d: alpha vs changed (alpha + beta)
plot_alpha_changed <- ggplot(
  data = alpha_beta_for_plotting %>% filter(bootstrapped == TRUE),
  aes(x = alpha_coefficients, y = changed)
) +
  geom_point(alpha = transparency, color = "grey") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(
    limits = c(
      quantile(
        (alpha_beta_for_plotting %>%
          filter(bootstrapped == TRUE))$alpha_coefficients,
        perc / 2
      ),
      quantile(
        (alpha_beta_for_plotting %>%
          filter(bootstrapped == TRUE))$alpha_coefficients,
        1 - (perc / 2)
      )
    )
  ) +
  scale_y_continuous(
    limits = c(
      quantile(
        (alpha_beta_for_plotting %>% filter(bootstrapped == TRUE))$changed,
        perc / 2
      ),
      quantile(
        (alpha_beta_for_plotting %>% filter(bootstrapped == TRUE))$changed,
        1 - (perc / 2)
      )
    )
  ) +
  geom_point(
    data = alpha_beta_for_plotting %>% filter(bootstrapped == FALSE),
    aes(x = alpha_coefficients, y = changed),
    color = "black",
    shape = 19
  ) +
  geom_smooth(
    data = alpha_beta_for_plotting %>% filter(bootstrapped == TRUE),
    aes(x = alpha_coefficients, y = changed),
    method = "glm",
    color = color_smooth,
    linetype = "dashed",
    se = FALSE
  ) +
  geom_smooth(
    data = alpha_beta_for_plotting %>% filter(bootstrapped == FALSE),
    aes(x = alpha_coefficients, y = changed),
    method = "glm",
    color = color_smooth,
    se = FALSE
  ) +
  xlab(expression("Plant-plant interactions (" * alpha * ")")) +
  ylab(expression("Plant-plant + higher-order (" * alpha + beta * ")")) +
  theme_classic()


# arrange panels c and d
arranged_hois <- ggpubr::ggarrange(
  plot_alpha_beta,
  plot_alpha_changed,
  ncol = 1,
  nrow = 2,
  align = "v",
  labels = c("c", "d"),
  vjust = 0.8
)


# final assembly ---------------------------------------------------------

# arrange everything into the final figure
arranged_all <- ggpubr::ggarrange(
  arranged_igr,
  NULL,
  arranged_hois,
  ncol = 3,
  nrow = 1,
  widths = c(0.5 - (0.075 / 2), 0.075, 0.5 - (0.075 / 2)),
  align = "hv"
)

# save the arranged plot
ggsave(
  arranged_all,
  file = "results/figures/fig-assess-complexity.jpeg",
  device = "jpeg",
  dpi = 320,
  height = 6.5,
  width = 7.75
)
