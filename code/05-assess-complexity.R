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
  boot_id = 0,
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

# 3. combine for scatter plot (panel a & b)
igr_gamma_for_plotting <- bind_rows(
  emp_igr_gamma %>% select(-spp),
  boot_igr_gamma
)

# create the net variable for panel b
igr_gamma_for_plotting <- igr_gamma_for_plotting %>%
  mutate(net_igr = igr + grasshoppers)


# data preparation: alpha and beta ---------------------------------------

# 1. extract empirical alpha and beta
emp_alpha_beta <- tibble(
  bootstrapped = FALSE,
  boot_id = 0,
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
  select(boot_id, alpha_coefficients, beta_coefficients) %>%
  unnest(c(alpha_coefficients, beta_coefficients)) %>%
  mutate(bootstrapped = TRUE)

# 3. combine for panels c and d
alpha_beta_for_plotting <- bind_rows(emp_alpha_beta, boot_alpha_beta) %>%
  mutate(changed = alpha_coefficients + beta_coefficients)


# stats function ---------------------------------------------------------

# function to calculate empirical R2 and bootstrap CI
get_r2_ci <- function(df_all, x_col, y_col) {
  # empirical model
  emp_df <- df_all %>% filter(bootstrapped == FALSE)
  formula_obj <- as.formula(paste(y_col, "~", x_col))
  m_emp <- lm(formula_obj, data = emp_df)
  # Updated to 3 decimal places
  emp_r2 <- round(summary(m_emp)$r.squared, 3)

  # bootstrap models
  boot_r2s <- df_all %>%
    filter(bootstrapped == TRUE) %>%
    group_by(boot_id) %>%
    summarise(
      r2 = summary(lm(formula_obj, data = pick(everything())))$r.squared,
      .groups = "drop"
    ) %>%
    pull(r2)

  # calculate quantiles
  if (length(unique(boot_r2s)) > 1) {
    ci <- quantile(boot_r2s, c(0.025, 0.975), na.rm = TRUE)
    # Updated to 3 decimal places
    ci_lab <- paste0("[", round(ci[1], 3), ", ", round(ci[2], 3), "]")
  } else {
    ci_lab <- "[NA, NA]"
  }

  # format label for ggplot parse
  label <- paste0("R^2 == ", emp_r2, "~'", ci_lab, "'")
  return(label)
}

# plotting limits setup --------------------------------------------------

# Panel A Y-axis limit
lims_a_y <- quantile(
  igr_gamma_for_plotting$grasshoppers,
  c(perc / 2, 1 - (perc / 2)),
  na.rm = TRUE
)

# Shared Panel A & B X-axis / Y-axis limits (forces vertical alignment)
lims_b <- c(
  min(
    quantile(igr_gamma_for_plotting$igr, perc / 2, na.rm = T),
    quantile(igr_gamma_for_plotting$net_igr, perc / 2, na.rm = T)
  ),
  max(
    quantile(igr_gamma_for_plotting$igr, 1 - perc / 2, na.rm = T),
    quantile(igr_gamma_for_plotting$net_igr, 1 - perc / 2, na.rm = T)
  )
)

# Panel C Y-axis limit
lims_c_y <- quantile(
  alpha_beta_for_plotting$beta_coefficients,
  c(perc / 2, 1 - (perc / 2)),
  na.rm = TRUE
)

# Shared Panel C & D X-axis / Y-axis limits (forces vertical alignment)
lims_d <- c(
  min(
    quantile(alpha_beta_for_plotting$alpha_coefficients, perc / 2, na.rm = T),
    quantile(alpha_beta_for_plotting$changed, perc / 2, na.rm = T)
  ),
  max(
    quantile(
      alpha_beta_for_plotting$alpha_coefficients,
      1 - perc / 2,
      na.rm = T
    ),
    quantile(alpha_beta_for_plotting$changed, 1 - perc / 2, na.rm = T)
  )
)


# plotting ---------------------------------------------------------------

# panel a: igr vs gamma
stats_a <- get_r2_ci(igr_gamma_for_plotting, "igr", "grasshoppers")

plot_r_gamma <- ggplot(
  data = igr_gamma_for_plotting %>% filter(bootstrapped == TRUE),
  aes(x = igr, y = grasshoppers)
) +
  geom_point(alpha = transparency, color = "grey") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_x_continuous(limits = lims_b) + # NOW SHARES LIMITS WITH PANEL B
  scale_y_continuous(limits = lims_a_y) +
  geom_point(
    data = igr_gamma_for_plotting %>% filter(bootstrapped == FALSE),
    aes(x = igr, y = grasshoppers),
    color = "black",
    shape = 19
  ) +
  geom_smooth(
    data = igr_gamma_for_plotting %>% filter(bootstrapped == FALSE),
    aes(x = igr, y = grasshoppers),
    method = "glm",
    color = color_smooth,
    se = FALSE
  ) +
  ggtitle(parse(text = stats_a)) +
  xlab(expression(atop("", "Intrinsic growth rates (r)"))) +
  ylab(expression(atop("Direct herbivory", "effects (" * gamma * ")"))) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 11, face = "plain", hjust = 0.5), # CENTERED TITLE
    axis.title.x = element_text(color = NA)
  )


# panel b: net effect on growth (igr vs igr + gamma)
stats_b <- get_r2_ci(igr_gamma_for_plotting, "igr", "net_igr")

plot_r_changed <- ggplot(
  data = igr_gamma_for_plotting %>% filter(bootstrapped == TRUE),
  aes(x = igr, y = net_igr)
) +
  geom_point(alpha = transparency, color = "grey") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_x_continuous(limits = lims_b) +
  scale_y_continuous(limits = lims_b) +
  geom_point(
    data = igr_gamma_for_plotting %>% filter(bootstrapped == FALSE),
    aes(x = igr, y = net_igr),
    color = "black",
    shape = 19
  ) +
  geom_smooth(
    data = igr_gamma_for_plotting %>% filter(bootstrapped == FALSE),
    aes(x = igr, y = net_igr),
    method = "glm",
    color = color_smooth,
    se = FALSE
  ) +
  ggtitle(parse(text = stats_b)) +
  xlab(expression(atop("", "Intrinsic growth rates (r)"))) +
  ylab(expression(atop(
    "Intrinsic growth rates",
    "with herbivory (r + " * gamma * ")"
  ))) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 11, face = "plain", hjust = 0.5) # CENTERED TITLE
  )


# panel c: alpha vs beta
stats_c <- get_r2_ci(
  alpha_beta_for_plotting,
  "alpha_coefficients",
  "beta_coefficients"
)

plot_alpha_beta <- ggplot(
  data = alpha_beta_for_plotting %>% filter(bootstrapped == TRUE),
  aes(x = alpha_coefficients, y = beta_coefficients)
) +
  geom_point(alpha = transparency, color = "grey") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(limits = lims_d) + # NOW SHARES LIMITS WITH PANEL D
  scale_y_continuous(limits = lims_c_y) +
  geom_point(
    data = alpha_beta_for_plotting %>% filter(bootstrapped == FALSE),
    aes(x = alpha_coefficients, y = beta_coefficients),
    color = "black",
    shape = 19
  ) +
  geom_smooth(
    data = alpha_beta_for_plotting %>% filter(bootstrapped == FALSE),
    aes(x = alpha_coefficients, y = beta_coefficients),
    method = "glm",
    color = color_smooth,
    se = FALSE
  ) +
  ggtitle(parse(text = stats_c)) +
  xlab(expression(atop("", "Plant-plant interactions (" * alpha * ")"))) +
  ylab(expression(atop("Higher-order interactions", "(HOIs; " * beta * ")"))) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 11, face = "plain", hjust = 0.5), # CENTERED TITLE
    axis.title.x = element_text(color = NA)
  )


# panel d: alpha vs changed (alpha + beta)
stats_d <- get_r2_ci(alpha_beta_for_plotting, "alpha_coefficients", "changed")

plot_alpha_changed <- ggplot(
  data = alpha_beta_for_plotting %>% filter(bootstrapped == TRUE),
  aes(x = alpha_coefficients, y = changed)
) +
  geom_point(alpha = transparency, color = "grey") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_x_continuous(limits = lims_d) +
  scale_y_continuous(limits = lims_d) +
  geom_point(
    data = alpha_beta_for_plotting %>% filter(bootstrapped == FALSE),
    aes(x = alpha_coefficients, y = changed),
    color = "black",
    shape = 19
  ) +
  geom_smooth(
    data = alpha_beta_for_plotting %>% filter(bootstrapped == FALSE),
    aes(x = alpha_coefficients, y = changed),
    method = "glm",
    color = color_smooth,
    se = FALSE
  ) +
  ggtitle(parse(text = stats_d)) +
  xlab(expression(atop("", "Plant-plant interactions (" * alpha * ")"))) +
  ylab(expression(atop(
    "Plant-plant interactions",
    "with HOIs (" * alpha + beta * ")"
  ))) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 11, face = "plain", hjust = 0.5) # CENTERED TITLE
  )


# final assembly ---------------------------------------------------------

# arrange everything into the final figure
arranged_all <- ggpubr::ggarrange(
  plot_r_gamma,
  NULL,
  plot_alpha_beta,
  plot_r_changed,
  NULL,
  plot_alpha_changed,
  ncol = 3,
  nrow = 2,
  align = "hv",
  widths = c(1, 0.15, 1),
  labels = c("a", "", "c", "b", "", "d"),
  hjust = -0.75,
  font.label = list(size = 14, face = "bold")
)

# save the arranged plot
ggsave(
  arranged_all,
  file = "results/figures/fig-3.jpeg",
  device = "jpeg",
  dpi = 320,
  height = 6.25,
  width = 7.25
)
