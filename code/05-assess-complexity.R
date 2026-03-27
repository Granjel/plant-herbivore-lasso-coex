# visual panels with linear regressions to assess the role of herbivory and HOIs

# setup ------------------------------------------------------------------

# load setup
source("code/01-setup.R")

# load empirical and bootstrapped matrices
load("data/processed/empirical/empirical-matrices.RData")
load("data/processed/bootstrapped/bootstrapped-matrices.RData")

# define ci to crop plot and regression lines color
perc <- 0.01
color_smooth <- c("darkgreen", "deepskyblue4")[1]
transparency <- 0.15


# data preparation: igr and gamma ----------------------------------------

# extract empirical igr and gamma
emp_igr_gamma <- tibble(
  spp_code = species,
  bootstrapped = FALSE,
  boot_id = 0,
  igr = igr,
  grasshoppers = rowSums(gamma)
)

# extract bootstrapped igr and gamma
boot_igr_gamma <- boot_matrices %>%
  dplyr::select(boot_id, igr, gamma) %>%
  mutate(
    spp_code = list(species),
    grasshoppers = map(gamma, rowSums)
  ) %>%
  dplyr::select(-gamma) %>%
  unnest(c(spp_code, igr, grasshoppers)) %>%
  mutate(bootstrapped = TRUE)

# combine for scatter plot (panel a & b)
igr_gamma_for_plotting <- bind_rows(emp_igr_gamma, boot_igr_gamma) %>%
  mutate(net_igr = igr + grasshoppers)


# data preparation: alpha and beta ---------------------------------------

# extract empirical alpha and beta
emp_alpha_beta <- tibble(
  bootstrapped = FALSE,
  boot_id = 0,
  alpha_coefficients = as.vector(alpha),
  beta_coefficients = as.vector(Reduce("+", beta_gp)) +
    as.vector(Reduce("+", beta_pp))
)

# extract bootstrapped alpha and beta
boot_alpha_beta <- boot_matrices %>%
  mutate(
    alpha_coefficients = map(alpha, as.vector),
    beta_coefficients = map2(
      beta_pp,
      beta_gp,
      ~ as.vector(Reduce("+", .x)) + as.vector(Reduce("+", .y))
    )
  ) %>%
  dplyr::select(boot_id, alpha_coefficients, beta_coefficients) %>%
  unnest(c(alpha_coefficients, beta_coefficients)) %>%
  mutate(bootstrapped = TRUE)

# combine for panels c and d
alpha_beta_for_plotting <- bind_rows(emp_alpha_beta, boot_alpha_beta) %>%
  mutate(changed = alpha_coefficients + beta_coefficients)


# stats function ---------------------------------------------------------

# function to calculate empirical R2 and bootstrap CI
get_r2_ci <- function(df_all, x_col, y_col) {
  # empirical model
  emp_df <- df_all %>% filter(bootstrapped == FALSE)
  formula_obj <- as.formula(paste(y_col, "~", x_col))
  m_emp <- lm(formula_obj, data = emp_df)
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
    ci <- quantile(boot_r2s, c(0.005, 0.995), na.rm = TRUE)
    ci_lab <- paste0("[", round(ci[1], 3), ", ", round(ci[2], 3), "]")
  } else {
    ci_lab <- "[NA, NA]"
  }

  # format label for ggplot parse
  label <- paste0("R^2 == ", emp_r2, "~'", ci_lab, "'")
  return(label)
}

# plotting limits setup --------------------------------------------------

# fn to keep all empirical data, plus the middle 98% of bootstrap data
get_smart_limits <- function(emp_vals, boot_vals, p = 0.01) {
  lower <- min(
    min(emp_vals, na.rm = TRUE),
    quantile(boot_vals, p, na.rm = TRUE)
  )
  upper <- max(
    max(emp_vals, na.rm = TRUE),
    quantile(boot_vals, 1 - p, na.rm = TRUE)
  )
  return(c(lower, upper))
}

# isolate data types
emp_ab <- igr_gamma_for_plotting %>% filter(bootstrapped == FALSE)
boot_ab <- igr_gamma_for_plotting %>% filter(bootstrapped == TRUE)

emp_cd <- alpha_beta_for_plotting %>% filter(bootstrapped == FALSE)
boot_cd <- alpha_beta_for_plotting %>% filter(bootstrapped == TRUE)

# panel a y-axis limit
lims_a_y <- get_smart_limits(emp_ab$grasshoppers, boot_ab$grasshoppers, perc)

# shared panel a & b limits
lims_b <- get_smart_limits(
  c(emp_ab$igr, emp_ab$net_igr),
  c(boot_ab$igr, boot_ab$net_igr),
  perc
)

# panel c y-axis limit
lims_c_y <- get_smart_limits(
  emp_cd$beta_coefficients,
  boot_cd$beta_coefficients,
  perc
)

# shared panel c & d limits
lims_d <- get_smart_limits(
  c(emp_cd$alpha_coefficients, emp_cd$changed),
  c(boot_cd$alpha_coefficients, boot_cd$changed),
  perc
)


# plotting ---------------------------------------------------------------

# panel a: igr vs gamma
stats_a <- get_r2_ci(igr_gamma_for_plotting, "igr", "grasshoppers")
plot_a <- ggplot(
  data = igr_gamma_for_plotting %>% filter(bootstrapped == TRUE),
  aes(x = igr, y = grasshoppers)
) +
  geom_point(alpha = transparency, color = "grey") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_x_continuous(limits = lims_b) +
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
    plot.title = element_text(size = 11, face = "plain", hjust = 0.5),
    axis.title.x = element_text(color = NA)
  )

# panel b: net effect on growth
stats_b <- get_r2_ci(igr_gamma_for_plotting, "igr", "net_igr")
plot_b <- ggplot(
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
  theme(plot.title = element_text(size = 11, face = "plain", hjust = 0.5))

# panel c: alpha vs beta
stats_c <- get_r2_ci(
  alpha_beta_for_plotting,
  "alpha_coefficients",
  "beta_coefficients"
)
plot_c <- ggplot(
  data = alpha_beta_for_plotting %>% filter(bootstrapped == TRUE),
  aes(x = alpha_coefficients, y = beta_coefficients)
) +
  geom_point(alpha = transparency, color = "grey") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(limits = lims_d) +
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
    plot.title = element_text(size = 11, face = "plain", hjust = 0.5),
    axis.title.x = element_text(color = NA)
  )

# panel d: alpha vs changed
stats_d <- get_r2_ci(alpha_beta_for_plotting, "alpha_coefficients", "changed")
plot_d <- ggplot(
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
  theme(plot.title = element_text(size = 11, face = "plain", hjust = 0.5))


# final assembly using patchwork -----------------------------------------

# patchwork formula to automatically handle layout and alignment
arranged_all <- (plot_a | plot_c) /
  (plot_b | plot_d) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 14, face = "bold"))

# save the arranged plot
ggsave(
  "results/figures/fig-complexity.jpeg",
  plot = arranged_all,
  device = "jpeg",
  dpi = dpi,
  height = 6.25,
  width = 7.25
)
