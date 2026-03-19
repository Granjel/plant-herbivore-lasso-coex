# set up specific for lasso (both empirical data and bootstraps)

# general setup ----------------------------------------------------------

# libraries, data, etc.
source("code/01-setup.R")

# create directory for empirical results if it doesn't exist
if (!dir.exists("data/processed/empirical")) {
  dir.create("data/processed/empirical")
}

# define columns for plants and grasshoppers
start_col_plants <- 8
end_col_plants <- 43 # plants start in col 8 and end in col 43
start_col_grasshoppers <- 47
end_col_grasshoppers <- 52 # cols 47 to 52 for grasshoppers

# define plant and grasshopper names
species <- colnames(data[start_col_plants:end_col_plants])
grasshoppers <- colnames(data[start_col_grasshoppers:end_col_grasshoppers])

# define the lambda values for lasso (default is 50)
max_lambda <- 50

# define subset for lasso
species_subset <- c(
  start_col_plants:end_col_plants,
  start_col_grasshoppers:end_col_grasshoppers
)

# define column count for lasso
n_col <- length(species_subset)
