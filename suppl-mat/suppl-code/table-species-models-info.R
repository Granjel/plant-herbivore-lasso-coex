# table with species codes, cover, individuals, and lambdas used

# setup ------------------------------------------------------------------

# load setup and libraries
source("code/01-setup.R")

# creating a df with individuals, cover, and species
data <- data %>%
  group_by(Focal) %>%
  summarise(
    indiv = n(),
    cover = sum(Cover)
  ) %>%
  rename(code = Focal) %>%
  arrange(code)

# calculate relative cover and individuals
data$indiv <- round(data$indiv / sum(data$indiv) * 100, 2)
data$cover <- round(data$cover / sum(data$cover) * 100, 2)

# add lambdas
load("data/processed/empirical/empirical-lambdas.RData")
data$lambdas <- round(lambda_table$lambda_value, 4)

# arrange
data <- data %>%
  mutate(abbrev = species_abbrev, full = species_full) %>%
  dplyr::select(code, abbrev, full, cover, indiv, lambdas) %>%
  arrange(-cover)

# save
write.table(
  data,
  file = "suppl-mat/suppl-tables/table-species-models-info.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
