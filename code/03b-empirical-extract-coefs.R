# extract empirical coefficients from glinternet and build matrices

# setup ------------------------------------------------------------------

# load lasso parameters
source("code/02-lasso-parameters.R")

# load lasso results
load("data/processed/empirical/empirical-lasso-models.RData")

# info
cat("Extracting coefficients and building interaction matrices...\n")


# extract lambdas from glinternet objects --------------------------------

# vectors to save selected lambdas
lambda_idx <- NULL # index along lambda path (what you currently store)
lambda_val <- NULL # actual lambda value

# select lambda values to ensure min error and intraspecific terms
for (i in 1:length(gli_models)) {
  # lambda selected for minimum error (index)
  n_i <- which(gli_models[[i]]$lambdaHat == gli_models[[i]]$lambda)

  # select lambda that also includes the intra
  for (j in n_i:max_lambda) {
    if (i %in% coef(gli_models[[i]]$glinternetFit)[[j]]$mainEffects$cont) {
      # store index
      lambda_idx <- c(lambda_idx, j)

      # store actual lambda value
      lambda_val <- c(lambda_val, gli_models[[i]]$lambda[j])

      break # stop j when lambda is found
    }

    if (j == max_lambda) {
      # fallback: store max index
      lambda_idx <- c(lambda_idx, max_lambda)

      # fallback: store corresponding lambda value
      lambda_val <- c(lambda_val, gli_models[[i]]$lambda[max_lambda])
    }
  } # end j
} # end i

# any lambda over the max?
if (isTRUE(length(which(lambda_idx >= max_lambda)) != 0)) {
  cat(
    "There is an issue: the max number of lambdas was not sufficient!\n"
  )
} else {
  cat("All good with the lambdas used!\n")
}

# combine into a table
lambda_table <- data.frame(
  species = species,
  lambda_index = lambda_idx,
  lambda_value = lambda_val
)

# save lambdas used for coefficient extraction
save(lambda_table, file = "data/processed/empirical/empirical-lambdas.RData")


# extract intrinsic growth rates, fixed, and interactive effects ---------

# objects to save growth rates, fixed effects, and HOIs
igr <- NULL # growth rates
fixed <- list() # intra & interspecific interactions
inter <- list() # HOIs

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
  file = "data/processed/empirical/empirical-coefficients.RData"
)
