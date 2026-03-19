# build the matrices for the empirical data

# setup ------------------------------------------------------------------

# load general setup and empirical setup
source("code/02-empirical-models.R")

# load lasso results
load("data/processed/empirical/empirical-coefficients-unproccessed.RData")

# info
cat("Building the interacion matrices with the empirical coefficients...\n")


# build alpha (community matrix) and gamma (herbivory matrix) ------------

# object to save both together first
alpha_gamma <- matrix(0, length(species), n_col)

# create a matrix for plant-plant + grasshopper-plant together
for (i in 1:length(species)) {
  for (j in 1:nrow(fixed[[i]])) {
    alpha_gamma[i, fixed[[i]]$pos[j]] <- fixed[[i]]$coef[j]
  }
}

# easier to record the position of plants and grasshoppers
plant_pos <- 1:length(species)
grasshopper_pos <- (length(species) + 1):n_col

# separate alpha (plant-plant pairwise interaction matrix)
alpha <- alpha_gamma[, plant_pos]
colnames(alpha) <- species
rownames(alpha) <- species

# separate gamma (grasshopper-plant pairwise interaction matrix)
gamma <- alpha_gamma[, grasshopper_pos]
colnames(gamma) <- grasshoppers
rownames(gamma) <- species


# build beta (HOIs) matrices ---------------------------------------------

# lists to save the different HOIs
inter_pp <- list() #plants on plant-plant
inter_gp <- list() #grasshoppers on plant-plant

# loop to separate the HOIs
for (i in 1:length(inter)) {
  if (nrow(inter[[i]]) > 0) {
    # plants on plant-plant
    inter_pp[[i]] <- inter[[i]] %>%
      filter((spp1 %in% plant_pos) & (spp2 %in% plant_pos))
    # grasshoppers on plant-plant
    inter_gp[[i]] <- inter[[i]] %>%
      filter((spp2 %in% grasshopper_pos) & (!spp1 %in% grasshopper_pos))
  } else {
    inter_pp[[i]] <- tibble()
    inter_gp[[i]] <- tibble()
  }
} #end i

# lists to build the actual matrices
beta_pp <- lapply(
  1:length(species),
  matrix,
  data = 0,
  nrow = length(species),
  ncol = length(species)
) #plants on plant-plant
beta_gp <- lapply(
  1:length(grasshoppers),
  matrix,
  data = 0,
  nrow = length(species),
  ncol = length(species)
) #grasshoppers on plant-plant

# loop to fill the matrices
for (i in 1:length(inter)) {
  # beta_pp
  for (j in 1:nrow(inter_pp[[i]])) {
    if (nrow(inter_pp[[i]]) > 0) {
      beta_pp[[inter_pp[[i]]$spp1[j]]][i, inter_pp[[i]]$spp2[j]] <- inter_pp[[
        i
      ]]$coef[j] /
        2
      beta_pp[[inter_pp[[i]]$spp2[j]]][i, inter_pp[[i]]$spp1[j]] <- inter_pp[[
        i
      ]]$coef[j] /
        2
    }
  }
  # beta_gp
  for (j in 1:nrow(inter_gp[[i]])) {
    if (nrow(inter_gp[[i]]) > 0) {
      beta_gp[[(inter_gp[[i]]$spp2[j] - length(species))]][
        i,
        inter_gp[[i]]$spp1[j]
      ] <- inter_gp[[i]]$coef[j]
    }
  }
}

# name beta_pp matrices
for (i in 1:length(beta_pp)) {
  rownames(beta_pp[[i]]) <- species
  colnames(beta_pp[[i]]) <- species
}

# name beta_gp matrices
for (i in 1:length(beta_gp)) {
  rownames(beta_gp[[i]]) <- species
  colnames(beta_gp[[i]]) <- species
}


# add plant names to the intrisic growth rate vector ---------------------

# name igr
names(igr) <- species


# save all matrices together ---------------------------------------------

# as an .RData file
save(
  igr, # r, or intrinsic growth rates of plants
  alpha, # alpha, or plant-plant pairwise interaction matrix
  gamma, # gamma, or grasshopper-plant pairwise interaction matrix
  beta_pp, # beta_pp, or plant-plant-plant HOIs matrix
  beta_gp, # beta_gp, or grasshopper-plant-plant HOIs matrix
  species, # keep plant species names
  grasshoppers, # keep grasshopper species names
  file = "data/processed/empirical/empirical-matrices.RData"
)
