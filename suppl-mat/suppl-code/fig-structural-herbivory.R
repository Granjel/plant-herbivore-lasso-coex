# evaluate the effect of herbivory on structural fitness differences (SFD)

# setup ------------------------------------------------------------------

# load setup
source("code/01-setup.R")
library(patchwork)

# load empirical matrices
load("data/processed/empirical/empirical-matrices.RData")

# prepare intrinsic growth rates (with and without herbivory)
names(igr) <- rownames(alpha)

# net growth rate = intrinsic + direct grasshopper effects
igr_herb <- igr + rowSums(gamma)
names(igr_herb) <- rownames(alpha)


# structural functions ---------------------------------------------------

# calculate the centroid of the feasibility domain
r_centroid <- function(alpha) {
  n <- nrow(alpha)
  D <- diag(1 / sqrt(diag(t(alpha) %*% alpha)))
  alpha_n <- alpha %*% D
  r_c <- rowSums(alpha_n) / n
  r_c <- t(t(r_c))
  return(r_c)
}

# structural fitness differences (SFD)
theta <- function(alpha, r) {
  r_c <- r_centroid(alpha)
  out <- acos(sum(r_c * r) / (sqrt(sum(r^2)) * sqrt(sum(r_c^2)))) * 180 / pi
  return(out)
}

# feasibility
test_feasibility <- function(alpha, r) {
  out <- prod(solve(alpha, r) > 0)
  return(out)
}

# compute wrapper specifically for herbivory comparison
eval_herbivory_sfd <- function(alpha, r_base, r_herb, n) {
  combos <- t(combn(rownames(alpha), n))

  res <- data.frame(
    combos = apply(combos, 1, paste, collapse = "_"),
    sfd_base = NA_real_,
    sfd_herb = NA_real_,
    feas_base = NA_integer_,
    feas_herb = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(combos)) {
    alpha2 <- as.matrix(alpha[combos[i, ], combos[i, ]])
    r_base2 <- as.vector(r_base[combos[i, ]])
    r_herb2 <- as.vector(r_herb[combos[i, ]])

    if (!all(alpha2 == 0)) {
      res$sfd_base[i] <- tryCatch(theta(alpha2, r_base2), error = function(e) {
        NA
      })
      res$sfd_herb[i] <- tryCatch(theta(alpha2, r_herb2), error = function(e) {
        NA
      })
      res$feas_base[i] <- tryCatch(
        test_feasibility(alpha2, r_base2),
        error = function(e) NA
      )
      res$feas_herb[i] <- tryCatch(
        test_feasibility(alpha2, r_herb2),
        error = function(e) NA
      )
    }
  }

  # calculate absolute shift here to keep plotting clean
  res <- res[complete.cases(res), ]
  res$sfd_shift <- res$sfd_herb - res$sfd_base
  return(res)
}


# computation ------------------------------------------------------------

cat("Computing SFD for 3-species modules...\n")
res_3 <- eval_herbivory_sfd(
  alpha = as.matrix(alpha),
  r_base = igr,
  r_herb = igr_herb,
  n = 3
)

cat("Computing SFD for 4-species modules (this may take a moment)...\n")
res_4 <- eval_herbivory_sfd(
  alpha = as.matrix(alpha),
  r_base = igr,
  r_herb = igr_herb,
  n = 4
)


# plotting ---------------------------------------------------------------

# helper function to keep plot aesthetics consistent
plot_sfd_hex <- function(df, legend_text) {
  ggplot(df, aes(x = sfd_base, y = sfd_shift)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_hex(bins = 50) +
    scale_fill_gradient(
      low = "grey90",
      high = "darkgreen",
      name = legend_text,
      guide = guide_colorbar(
        barwidth = unit(3, "cm"), # prevents text overlap
        barheight = unit(0.25, "cm") # <-- Makes the gradient rectangle slimmer!
      )
    ) +
    labs(
      x = "Structural fitness differences (SFD)",
      y = expression("Shift due to herbivory (" * Delta * "SFD)")
    ) +
    theme_classic() +
    theme(
      legend.position = "top",
      legend.title = element_text(lineheight = 1.2) # keeps the \n spacing clean
    )
}

# create the individual plots
plot_3_sp <- plot_sfd_hex(res_3, "3-species modules (count)") +
  theme(axis.title.x = element_blank()) # hide x-axis title for the top plot

plot_4_sp <- plot_sfd_hex(res_4, "4-species modules (count)")

# combine them vertically and make tags bold
final_combined_plot <- (plot_3_sp / plot_4_sp) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 14, face = "bold"))


# final export -----------------------------------------------------------

# save combined plot
ggsave(
  filename = "suppl-mat/suppl-figures/fig-sfd-herbivory.jpeg",
  plot = final_combined_plot,
  device = "jpeg",
  dpi = 640,
  height = 7.75,
  width = 4
)

# info
cat("Done! Combined vertical plot and data tables saved to suppl-mat/.\n")
