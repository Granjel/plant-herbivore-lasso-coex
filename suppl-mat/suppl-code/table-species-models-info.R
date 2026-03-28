# table with species codes, cover, individuals, and lambdas used

# setup ------------------------------------------------------------------

# load setup and libraries
source("code/01-setup.R")

# load full community data
full_community <- read.table(
  "data/raw/full-plant-community-cover.txt",
  header = TRUE,
  sep = "\t"
)

# species codes
species_codes <- data %>%
  pull(Focal) %>%
  unique

# get total and % cover per species in the full community
cover_per_species <- full_community %>%
  group_by(Spcode) %>%
  summarise(cover_absolute = sum(Cover), .groups = "drop") %>%
  mutate(cover_percentage = cover_absolute / sum(cover_absolute) * 100)

# get total and % individuals per species in the full community
individuals_per_species <- full_community %>%
  group_by(Spcode) %>%
  summarise(individuals_absolute = n(), .groups = "drop") %>%
  mutate(
    individuals_percentage = individuals_absolute /
      sum(individuals_absolute) *
      100
  )

# create a df with all species, percentage cover, and percentage individuals
df <- data.frame(
  code = cover_per_species$Spcode,
  cover = round(cover_per_species$cover_percentage, 2),
  indiv = round(individuals_per_species$individuals_percentage, 2)
)

# order by percentage cover and select the species that make up 98.75% of cover
df <- df %>%
  arrange(-cover) %>%
  mutate(cumulative_cover = cumsum(cover)) %>%
  filter(cumulative_cover <= 98.75) %>%
  dplyr::select(-cumulative_cover)

# add lambdas by left joining
load("data/processed/empirical/empirical-lambdas.RData")
df <- df %>%
  left_join(
    lambda_table %>%
      rename(code = species) %>%
      dplyr::select(code, lambda_value)
  )

# all types of species names
spp <- data.frame(
  code = species_codes
) %>%
  arrange(code) %>%
  mutate(
    abbrev = species_abbrev,
    full = species_full
  )

# join all, arrange and round
df <- spp %>%
  left_join(df, by = "code") %>%
  arrange(-cover) %>%
  mutate(lambda_value = round(lambda_value, 4))

# save
write.table(
  df,
  file = "suppl-mat/suppl-tables/table-species-models-info.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
