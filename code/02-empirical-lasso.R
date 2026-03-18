# empirical interaction coefficients with LV and lasso regularisation

# load setup -------------------------------------------------------------

# libraries, data, etc.
source("code/01-setup.R")

# define columns for plants and grasshoppers
start_col_plants <- 8
end_col_plants <- 43 #plants start in col 8 and end in col 43
start_col_grasshoppers <- 47
end_col_grasshoppers <- 52 #cols 47 to 52 for grasshoppers

# define plant and grasshopper names
species <- colnames(data[start_col_plants:end_col_plants])
grasshoppers <- colnames(data[start_col_grasshoppers:end_col_grasshoppers])

# define the lambda values for lasso (default is 50)
max_lambda <- 50

# create a list to save the lasso models from glinternet
gli_models <- list()

# define subset for lasso and column count
species_subset <- c(
  start_col_plants:end_col_plants,
  start_col_grasshoppers:end_col_grasshoppers
)
n_col <- length(c(
  start_col_plants:end_col_plants,
  start_col_grasshoppers:end_col_grasshoppers
))

## lasso: run only once!
# timer start
start <- Sys.time()

# lasso loop; takes some time
for (i in 1:length(species)) {
  gli_models[[i]] <- glinternet.cv(
    X = data[data$Focal == species[i], ][, species_subset],
    Y = data[data$Focal == species[i], ]$Cover,
    numLevels = rep(1, n_col),
    nLambda = max_lambda
  )
  cat(round(i / length(species) * 100, 1), "%\n")
} #end i

# timer end and duration
end <- Sys.time()
end - start

# save a .RData file with results after computing lasso
save(gli_models, file = "results/empirical/empirical-glinternet-models.RData")

# load lasso results
load("results/empirical/empirical-glinternet-models.RData")

# select lambda values to ensure min error and intraspecific terms; takes a little time
lambdas <- NULL #vector to save selected lambdas
for (i in 1:length(gli_models)) {
  # lambda selected for minimum error
  n_i <- which(gli_models[[i]]$lambdaHat == gli_models[[i]]$lambda)
  # select lambda that also includes the intra
  for (j in n_i:max_lambda) {
    if (i %in% coef(gli_models[[i]]$glinternetFit)[[j]]$mainEffects$cont) {
      lambdas <- c(lambdas, j)
      break #stop j when lambda is found
    }
    if (j == max_lambda) {
      lambdas <- c(lambdas, max_lambda)
    }
  } #end j
} #end i

# any lambda over the max?
if (isTRUE(length(which(lambdas >= 50)) != 0)) {
  cat(
    "MESSAGE: There is an issue: the max number of lambdas was not sufficient; check affected species!!!\n"
  )
} else {
  cat("MESSAGE: All good with the lambdas used!\n")
}

# save lambdas
save(lambdas, file = "results/empirical/empirical-lambdas.RData")

# objects to save growth rates, fixed effects, and HOIs
igr <- NULL #growth rates
fixed <- list() #intra & interspecific interactions
inter <- list() #HOIs

# loop to save such parameters
for (i in 1:length(gli_models)) {
  # select lambda to extract coefficients
  n_i <- lambdas[i]
  # extract corresponding coefficients
  coefs <- coef(gli_models[[i]]$glinternetFit)[[n_i]]
  # save IGR
  igr <- c(igr, gli_models[[i]]$betahat[[1]][1])
  # save fixed effects
  fixed[[i]] <- data.frame(
    "pos" = coefs$mainEffects$cont,
    "coef" = unlist(coefs$mainEffectsCoef$cont)
  )
  if (nrow(fixed[[i]]) > 0) {
    fixed[[i]] <- fixed[[i]][order(fixed[[i]]$pos), ]
  }
  # save HOIs
  inter[[i]] <- data.frame(
    "spp1" = coefs$interactions$contcont[, 1],
    "spp2" = coefs$interactions$contcont[, 2],
    "coef" = unlist(coefs$interactionsCoef$contcont)
  )
} #end i

# save a .RData file with extracted coefficients after computing lasso
save(
  species,
  grasshoppers,
  igr,
  fixed,
  inter,
  file = "results/empirical/empirical-coefficients.RData"
)

# create a matrix for plant-plant + grasshopper-plant together
alpha_gamma <- matrix(0, length(species), n_col)
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

# name igr
names(igr) <- species

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

# save all results
save(
  igr,
  alpha,
  gamma,
  beta_pp,
  beta_gp,
  species,
  grasshoppers,
  file = "results/empirical/empirical-matrices.RData"
)
