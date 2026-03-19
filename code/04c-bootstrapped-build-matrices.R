# build the matrices for the bootstrapped data

# setup ------------------------------------------------------------------

# load lasso parameters
source("code/02-lasso-parameters.R")

# load lasso results
load("data/processed/empirical/empirical-coefficients-unproccessed.RData")

# info
cat("Building the interacion matrices with the empirical coefficients...\n")
