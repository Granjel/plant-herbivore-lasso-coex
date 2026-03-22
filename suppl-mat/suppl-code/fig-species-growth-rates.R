# supplementary figure: forest plot of species-specific growth rates

# setup ------------------------------------------------------------------

# load empirical and bootstrapped matrices
# (assuming working directory is the project root and 01-setup.R has been run)
load("data/processed/empirical/empirical-matrices.RData")
load("data/processed/bootstrapped/bootstrapped-matrices.RData")

# full species names (mapped to the exact order of the original species vector)
species_full <- c(
  "Achillea millefolium",
  "Anthoxanthum odoratum",
  "Arrhenatherum elatius",
  "Bromus erectus",
  "Centaurea jacea",
  "Convolvulus arvensis",
  "Crepis sp.",
  "Dactylis glomerata",
  "Daucus carota",
  "Elytrigia repens",
  "Eryngium sp.",
  "Festuca arundinacea",
  "Festuca rubra",
  "Galium verum",
  "Geranium dissectum",
  "Geranium rotundifolium",
  "Leucanthemum vulgare",
  "Lolium perenne",
  "Lotus corniculatus",
  "Medicago arabica",
  "Ononis repens",
  "Picris echioides",
  "Picris hieracioides",
  "Plantago lanceolata",
  "Poa angustifolia",
  "Poa pratensis",
  "Poa trivialis",
  "Ranunculus acris",
  "Rumex acetosa",
  "Salvia pratensis",
  "Sonchus asper",
  "Taraxacum officinale",
  "Trifolium fragiferum",
  "Trifolium pratense",
  "Verbena officinalis",
  "Veronica persica"
)

# create plotmath-ready labels (italicize everything except 'sp.')
parse_labels <- paste0("italic('", species_full, "')")
parse_labels <- str_replace(parse_labels, " sp.'\\)", "')~plain('sp.')")
names(parse_labels) <- species_full

# data preparation -------------------------------------------------------

# extract empirical igr and gamma
emp_igr_gamma <- tibble(
  spp_code = species,
  spp = species_full,
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
  unnest(c(spp_code, igr, grasshoppers))

# calculate 95% CIs and format for forest plot
igr_forest_data <- boot_igr_gamma %>%
  group_by(spp_code) %>%
  summarise(
    igr_low = quantile(igr, 0.025),
    igr_high = quantile(igr, 0.975),
    igr_gamma_low = quantile(igr + grasshoppers, 0.025),
    igr_gamma_high = quantile(igr + grasshoppers, 0.975),
    .groups = "drop"
  ) %>%
  left_join(
    emp_igr_gamma %>%
      dplyr::select(spp_code, spp, emp_igr = igr, emp_g = grasshoppers),
    by = "spp_code"
  ) %>%
  uncount(2, .id = "condition") %>%
  mutate(
    # Create a factor with explicit levels so 'Without' gets dodged to the TOP
    added_gamma = factor(
      ifelse(condition == 1, "Without", "With"),
      levels = c("With", "Without")
    ),
    plot_igr = ifelse(added_gamma == "With", emp_igr + emp_g, emp_igr),
    plot_low = ifelse(added_gamma == "With", igr_gamma_low, igr_low),
    plot_high = ifelse(added_gamma == "With", igr_gamma_high, igr_high)
  )

# plotting ---------------------------------------------------------------

plot_forest <- ggplot(
  igr_forest_data,
  aes(
    y = reorder(spp, -plot_igr, FUN = median),
    x = plot_igr,
    shape = added_gamma,
    fill = added_gamma
  )
) +
  geom_pointrange(
    aes(xmin = plot_low, xmax = plot_high),
    position = position_dodge(width = 0.5),
    size = 0.3,
    color = "black",
    stroke = 0.5
  ) +
  scale_shape_manual(
    name = c("Intrinsic growth rates"),
    breaks = c("Without", "With"),
    values = c("Without" = 19, "With" = 21),
    labels = c(
      "Without" = "Without herbivory (r)",
      "With" = expression("With herbivory (" * r + gamma * ")")
    )
  ) +
  scale_fill_manual(
    name = c("Intrinsic growth rates"),
    breaks = c("Without", "With"),
    values = c("Without" = "black", "With" = "white"),
    labels = c(
      "Without" = "Without herbivory (r)",
      "With" = expression("With herbivory (" * r + gamma * ")")
    )
  ) +
  scale_y_discrete(labels = function(x) parse(text = parse_labels[x])) +
  xlab("Value") +
  ylab("Plant species") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = c(0.675, 0.875),
    legend.background = element_blank()
  )

# save the plot
ggsave(
  plot_forest,
  file = "suppl-mat/suppl-figures/fig-species-growth-rates.jpeg",
  device = "jpeg",
  dpi = 320,
  height = 7.0,
  width = 5.5
)
