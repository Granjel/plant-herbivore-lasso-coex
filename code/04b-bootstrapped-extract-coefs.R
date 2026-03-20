# extract bootstrapped coefficients and balance dropped iterations

# setup ------------------------------------------------------------------

# load lasso parameters
source("code/02-lasso-parameters.R")

# load bootstraps
load("data/processed/bootstrapped/bootstrapped-lasso-models.RData")

# rename object and save number of rows
number_of_rows <- nrow(bootstraps)

# define new columns
bootstraps$igr <- NA
bootstraps$lambda <- NA
bootstraps$delete <- FALSE
bootstraps$fixed <- lapply(
  1:number_of_rows,
  function(x) tibble(pos = numeric(), coef = numeric())
)
bootstraps$inter <- lapply(
  1:number_of_rows,
  function(x) tibble(spp1 = numeric(), spp2 = numeric(), coef = numeric())
)

# extract coefficients ---------------------------------------------------

# info
cat("Extracting coefficients...\n")

# loop to extract parameters and map to absolute positions
for (i in 1:number_of_rows) {
  # load model info
  model_info <- bootstraps$gli_models[[i]]

  # handle extinct focal species
  if (model_info$status == "extinct") {
    bootstraps$igr[i] <- 0
    next
  }

  # extract fit and index map
  fit <- model_info$fit
  index_map <- model_info$index_map

  # save the intrinsic growth rate
  bootstraps$igr[i] <- fit$betahat[[1]][1]

  # select lambda
  k_orig <- which(species == bootstraps$Focal[i])
  k_shrunken <- which(index_map == k_orig)

  n_i <- which(fit$lambdaHat == fit$lambda)

  for (j in n_i:max_lambda) {
    # from min error lambda until the intraspecific coefficient is found
    if (k_shrunken %in% coef(fit$glinternetFit)[[j]]$mainEffects$cont) {
      break # stop j when the desired lambda is found
    }
  } # end j

  # save lambda value
  bootstraps$lambda[i] <- j

  # identify those that didn't report an intraspecific value
  if (!k_shrunken %in% coef(fit$glinternetFit)[[j]]$mainEffects$cont) {
    bootstraps$delete[i] <- TRUE
  }

  # extract corresponding coefficients
  coefs <- coef(fit$glinternetFit)[[j]]

  # save fixed effects using the index map for absolute positioning
  bootstraps$fixed[[i]] <- tibble(
    pos = index_map[coefs$mainEffects$cont],
    coef = unlist(coefs$mainEffectsCoef$cont)
  )

  # arrange fixed effects by species number
  if (nrow(bootstraps$fixed[[i]]) > 0) {
    bootstraps$fixed[[i]] <- bootstraps$fixed[[i]][
      order(bootstraps$fixed[[i]]$pos),
    ]
  }

  # save HOIs using the index map for absolute positioning
  bootstraps$inter[[i]] <- tibble(
    spp1 = index_map[coefs$interactions$contcont[, 1]],
    spp2 = index_map[coefs$interactions$contcont[, 2]],
    coef = unlist(coefs$interactionsCoef$contcont)
  )
}

# alarm
beepr::beep(1)


# balance the dropped bootstraps ----------------------------------------

# info
cat("Balancing deleted iterations...\n")

# set seed for reproducible random deletions
set.seed(1610)

# count the number of TRUE values for each species in the delete column
true_counts <- bootstraps %>%
  group_by(Focal) %>%
  summarize(true_count = sum(delete == TRUE))

# find the maximum number of TRUE values across all species
max_true_count <- max(true_counts$true_count)

# for each species, randomly assign additional TRUE values to the delete column
bootstraps_balanced <- bootstraps %>%
  group_by(Focal) %>%
  mutate(
    # count current TRUE values for each species
    current_true_count = sum(delete == TRUE),
    # calculate how many more TRUE values are needed
    needed_true_count = max_true_count - current_true_count,
    # assign additional random deletions to balance datasets
    delete = ifelse(
      delete == TRUE,
      TRUE,
      ifelse(
        needed_true_count > 0 & sum(delete == FALSE) > 0,
        ifelse(
          row_number() %in%
            sample(
              which(delete == FALSE),
              min(needed_true_count, sum(delete == FALSE)),
              replace = FALSE
            ),
          TRUE,
          FALSE
        ),
        delete
      )
    )
  ) %>%
  ungroup()

# assign new random delete values and remove unnecessary objects
bootstraps$delete <- bootstraps_balanced$delete
rm(bootstraps_balanced, true_counts)

# delete those rows with delete == TRUE
bootstraps <- bootstraps %>% filter(delete == FALSE)

# info
cat(
  "Remaining successful bootstraps per species:",
  nrow(bootstraps) / length(species),
  "\n"
)


# structure for matrix building ------------------------------------------

# assign a row number within species, then order identically to the 'species' vector
bootstraps <- bootstraps %>%
  group_by(Focal) %>%
  mutate(boot_id = row_number()) %>%
  ungroup() %>%
  mutate(Focal = factor(Focal, levels = species)) %>%
  arrange(boot_id, Focal)

# new tibble to save all interaction matrices
boot_coefs <- bootstraps %>%
  group_by(boot_id) %>%
  summarise(
    Focal = list(Focal),
    fixed = list(fixed),
    inter = list(inter),
    igr = list(igr)
  ) %>%
  ungroup()

# clean up and save
rm(bootstraps)
save(
  boot_coefs,
  file = "data/processed/bootstrapped/bootstrapped-coefficients.RData"
)
