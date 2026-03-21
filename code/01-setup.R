# set up for analyses (packages, data, etc.)

# R libraries ------------------------------------------------------------

library(tidyverse) # for data manipulation and visualisation
library(glinternet) # for lasso regression with interactions
library(doSNOW) # for parallel processing with foreach and glinternet
library(foreach) # for parallel processing with doSNOW and glinternet
library(beepr) # for sound notification when done
library(patchwork) # for arranging ggplots


# create main directories ------------------------------------------------

# for processed or intermediate data
if (!dir.exists("data/processed")) {
  dir.create("data/processed")
}

# for results
if (!dir.exists("results")) {
  dir.create("results")
  dir.create("results/figures")
  dir.create("results/tables")
}

# for supplementary materials
if (!dir.exists("suppl-mat")) {
  dir.create("suppl-mat")
  dir.create("suppl-mat/suppl-code")
  dir.create("suppl-mat/suppl-figures")
  dir.create("suppl-mat/suppl-tables")
}


# theme for figures and colors -------------------------------------------

# overall theme for figures
theme_set(theme_bw(base_size = 14))


# empirical data ---------------------------------------------------------

# load plant-grasshopper data
data <- read.table("data/raw/empirical-dataset.txt", header = TRUE)
