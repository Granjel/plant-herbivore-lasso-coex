# build interaction matrices from bootstrapped coefficients

# setup ------------------------------------------------------------------

# load lasso parameters
source("code/02-lasso-parameters.R")

# load bootstrapped coefficients
load("data/processed/bootstrapped/bootstrapped-coefficients.RData")

# info
cat("Building the interaction matrices with the bootstrapped coefficients...\n")

# easier to record the position of plants and grasshoppers
plant_pos <- 1:length(species)
grasshopper_pos <- (length(species) + 1):n_col

# define new columns for matrices
boot_coefs$alpha <- lapply(1:nrow(boot_coefs), function(x) list())
boot_coefs$gamma <- lapply(1:nrow(boot_coefs), function(x) list())
boot_coefs$beta_pp <- lapply(1:nrow(boot_coefs), function(x) list())
boot_coefs$beta_gp <- lapply(1:nrow(boot_coefs), function(x) list())


# build alpha (community matrix) and gamma (herbivory matrix) ------------

# loop to build the matrices for each bootstrap
for (b in 1:nrow(boot_coefs)) {
  # object to save both together first
  alpha_gamma <- matrix(0, length(species), n_col)

  # create a matrix for plant-plant + grasshopper-plant together
  for (i in 1:length(species)) {
    if (nrow(boot_coefs$fixed[[b]][[i]]) > 0) {
      for (j in 1:nrow(boot_coefs$fixed[[b]][[i]])) {
        alpha_gamma[i, boot_coefs$fixed[[b]][[i]]$pos[j]] <- boot_coefs$fixed[[
          b
        ]][[
          i
        ]]$coef[j]
      }
    }
  }

  # separate alpha (plant-plant pairwise interaction matrix)
  boot_coefs$alpha[[b]] <- alpha_gamma[, plant_pos]
  colnames(boot_coefs$alpha[[b]]) <- species
  rownames(boot_coefs$alpha[[b]]) <- species

  # separate gamma (grasshopper-plant pairwise interaction matrix)
  boot_coefs$gamma[[b]] <- alpha_gamma[, grasshopper_pos]
  colnames(boot_coefs$gamma[[b]]) <- grasshoppers
  rownames(boot_coefs$gamma[[b]]) <- species

  # build beta (HOIs) matrices -------------------------------------------

  # lists to save the different HOIs
  inter_pp <- list() #plants on plant-plant
  inter_gp <- list() #grasshoppers on plant-plant

  # loop to separate the HOIs
  for (i in 1:length(boot_coefs$inter[[b]])) {
    if (nrow(boot_coefs$inter[[b]][[i]]) > 0) {
      # plants on plant-plant
      inter_pp[[i]] <- boot_coefs$inter[[b]][[i]] %>%
        filter((spp1 %in% plant_pos) & (spp2 %in% plant_pos))
      # grasshoppers on plant-plant
      inter_gp[[i]] <- boot_coefs$inter[[b]][[i]] %>%
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
  for (i in 1:length(boot_coefs$inter[[b]])) {
    # beta_pp
    if (nrow(inter_pp[[i]]) > 0) {
      for (j in 1:nrow(inter_pp[[i]])) {
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
    if (nrow(inter_gp[[i]]) > 0) {
      for (j in 1:nrow(inter_gp[[i]])) {
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

  # assign final beta matrices to the dataframe
  boot_coefs$beta_pp[[b]] <- beta_pp
  boot_coefs$beta_gp[[b]] <- beta_gp

  # add plant names to the intrisic growth rate vector -------------------

  names(boot_coefs$igr[[b]]) <- species
} # end b loop


# clean up and save ------------------------------------------------------

# leave only what's of interest
boot_matrices <- boot_coefs %>%
  dplyr::select(-Focal, -fixed, -inter)

# clean environment
rm(boot_coefs, alpha_gamma, beta_pp, beta_gp, inter_pp, inter_gp, i, j, b)

# save matrices
save(
  boot_matrices,
  file = "data/processed/bootstrapped/bootstrapped-matrices.RData"
)

# info
cat("Bootstrapped matrices successfully built and saved!\n")
