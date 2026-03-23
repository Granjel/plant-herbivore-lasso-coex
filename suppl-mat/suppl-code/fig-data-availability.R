# evaluate empirical data availability for interaction coefficients
# generates a 2x4 panel of specific focal ~ neighbor scatter plots

# setup ------------------------------------------------------------------

# load setup
source("code/01-setup.R")

# identify plant columns dynamically
start_col_plants <- which(colnames(data) == "ACHMIL")
end_col_plants <- which(colnames(data) == "VERPER")
species <- colnames(data)[start_col_plants:end_col_plants]

# clean and filter common tags across timepoints
data_clean <- data %>%
  filter(time %in% c(1, 3, 6)) %>%
  dplyr::select(
    time,
    block,
    treatment,
    datapoint,
    focal = Focal,
    cover = Cover,
    all_of(species)
  ) %>%
  mutate(tag = paste(block, treatment, datapoint, focal))

# identify tags that are present in all three time points
common_tags <- data_clean %>%
  distinct(time, tag) %>%
  group_by(time) %>%
  summarise(tags = list(tag)) %>%
  pull(tags) %>%
  Reduce(intersect, .)

# filter data to only include these common tags and arrange by time and tag
data_clean <- data_clean %>%
  filter(tag %in% common_tags) %>%
  arrange(time, tag)

# compute population change between time points
data1 <- data_clean %>%
  filter(time == 1) %>%
  mutate(time = 1, cover = cover + 1) %>%
  arrange(tag)
data2 <- data_clean %>%
  filter(time == 3) %>%
  mutate(time = 2, cover = cover + 1) %>%
  arrange(tag)
data3 <- data_clean %>%
  filter(time == 6) %>%
  mutate(time = 3, cover = cover + 1) %>%
  arrange(tag)

# calculate change in cover from time 1 to time 3, and from time 3 to time 6
data1 <- data1 %>% mutate(change = data2$cover / cover)
data2 <- data2 %>% mutate(change = data3$cover / cover)

# combine the two datasets and filter out any rows where change is NA
data_change <- bind_rows(data1, data2) %>% filter(!is.na(change))


# targeted species pairs -------------------------------------------------

# explicitly define the 8 pairs shown in your original figure
target_pairs <- list(
  # intraspecific (row 1-2)
  list(focal = "BROERE", neigh = "BROERE"),
  list(focal = "TRIPRA", neigh = "TRIPRA"),
  list(focal = "FESARU", neigh = "FESARU"),
  list(focal = "TAROFF", neigh = "TAROFF"),

  # interspecific (row 3-4)
  list(focal = "BROERE", neigh = "ARRELA"),
  list(focal = "ARRELA", neigh = "LOTCOR"),
  list(focal = "PLALAN", neigh = "POATRI"),
  list(focal = "ERYNGE", neigh = "TRIPRA")
)


# generate plots ---------------------------------------------------------

# function to create a scatter plot with a fitted line for each focal-neighbor
plot_list <- lapply(target_pairs, function(p) {
  # filter data for this specific interaction
  dat <- data_change %>%
    filter(focal == p$focal) %>%
    filter(.data[[p$neigh]] > 0)

  # explicitly ensure we have data before passing to lm()
  if (nrow(dat) > 2) {
    # safe formula build
    f <- as.formula(paste("change ~", p$neigh))
    mod <- lm(f, data = dat)
    dat$pred <- predict(mod)

    ggplot(dat, aes(x = .data[[p$neigh]], y = change)) +
      geom_point(
        position = position_jitter(width = 0.1, height = 0.1),
        color = "black",
        alpha = 0.75
      ) +
      geom_line(aes(y = pred), color = "forestgreen", linewidth = 1) +
      labs(
        x = paste(p$neigh, "cover"),
        y = paste(p$focal, "\ncover change")
      ) +
      theme_classic() +
      theme(
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)
      )
  } else {
    # fail-safe just in case
    ggplot() + theme_void() + ggtitle(paste("No data:", p$focal, "~", p$neigh))
  }
})


# assemble and export ----------------------------------------------------

# arrange the 8 plots into a 2x4 grid with labels
selected_all <- ggpubr::ggarrange(
  plotlist = plot_list,
  labels = "auto",
  align = "hv",
  ncol = 2,
  nrow = 4,
  font.label = list(size = 16, face = "bold")
)

# save the figure
ggsave(
  filename = "suppl-mat/suppl-figures/fig-data-availability.jpeg",
  plot = selected_all,
  device = "jpeg",
  width = 8,
  height = 11,
  dpi = 400,
  bg = "white"
)

# info
cat("Done! Targeted species comparisons successfully plotted.\n")
