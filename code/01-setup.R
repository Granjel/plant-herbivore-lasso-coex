# set up for analyses (packages, data, etc.)

# R libraries ------------------------------------------------------------

library(tidyverse) # for data manipulation and visualisation
library(glinternet) # for lasso regression with interactions
library(doSNOW) # for parallel processing with foreach and glinternet
library(foreach) # for parallel processing with doSNOW and glinternet
library(beepr) # for sound notification when done
library(patchwork) # for arranging ggplots
library(mvtnorm) # for multivariate normal distribution in str stab analysis
library(MASS) # for various statistical functions
library(EnvStats) # for environmental statistics
library(broom) # for tidying model outputs
library(igraph) # for network analysis and visualization
library(cowplot) # for arranging ggplots and saving figures
library(magick) # for image manipulation and saving figures

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
  dir.create("suppl-mat/packages")
}


# theme for figures and colors -------------------------------------------

# overall theme for figures
theme_set(theme_bw(base_size = 14))

# dpi for saving figures
dpi <- 640


# empirical data ---------------------------------------------------------

# load plant-grasshopper data
data <- read.table("data/raw/empirical-dataset.txt", header = TRUE)


# species names (abbreviation, full) -------------------------------------

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
  "En",
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

# grateful to cite packages used in the analyses
library(grateful)
cite_packages(
  out.format = "docx",
  out.dir = "suppl-mat/packages",
  citation.style = "oikos"
)
