# bootstrapped interaction coefficients with LV and lasso regularisation

# setup and compute bootstrapped data ------------------------------------

# load lasso parameters
source("code/02-lasso-parameters.R")

# times to resample (number of bootstraps)
t_boot <- 100

# nest the data per species
df <- data %>%
  group_by(Focal) %>%
  nest()

# function for resampling (bootstrap)
resample_the_data <- function(i, dat) {
  n_row <- nrow(dat)
  # Sample rows with replacement
  new_dat <- dat[sample(1:n_row, n_row, replace = TRUE), ]
  return(new_dat)
}

# set seed for reproducible resampling
set.seed(1610)

# resample the data within species
bootstraps <- df %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    resampled_data = list(lapply(1:t_boot, resample_the_data, dat = data))
  ) %>%
  unnest_longer(resampled_data)


# setup parallel computation ---------------------------------------------

# setup the cluster
n_cores <- parallel::detectCores(logical = TRUE) # - 2 # leave 2 cores free
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# assign a distinct random stream to each core for perfect reproducibility
parallel::clusterSetRNGStream(cl, iseed = 16010)

# setup the progress bar
total_tasks <- nrow(bootstraps)
pb <- txtProgressBar(max = total_tasks, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


# run bootstrapped lasso regularisation ----------------------------------

# info and timer start
cat("Running bootstrapped LASSO in parallel on", n_cores, "cores...\n")
cat("Total tasks to compute:", total_tasks, "\n")
start <- Sys.time()

# parallel loop
lasso_results <- foreach(
  i = 1:total_tasks,
  .packages = "glinternet",
  .options.snow = opts
) %dopar%
  {
    # extract the specific resampled dataframe for this iteration
    dat_i <- bootstraps$resampled_data[[i]]

    # check for ecological extinction of the focal species
    if (sum(dat_i$Cover) == 0) {
      return(list(status = "extinct"))
    }

    # identify present neighbours to avoid zero-variance columns
    X_full <- dat_i[, c(species, grasshoppers)]
    present_cols <- colSums(X_full) > 0
    X_present <- X_full[, present_cols]

    # run the model on the surviving community
    fit <- glinternet.cv(
      X = X_present,
      Y = dat_i$Cover,
      numLevels = rep(1, ncol(X_present)),
      nLambda = max_lambda
    )

    # return model and index map for proper coefficient placement
    return(list(
      status = "success",
      fit = fit,
      index_map = which(present_cols)
    ))
  }

# stop the cluster and close the progress bar
close(pb)
stopCluster(cl)

# timer end
end <- Sys.time()
print(end - start)

# sound notification when done
beepr::beep(1)


# save results -----------------------------------------------------------

# bind the results back to the tibble
bootstraps$gli_models <- lasso_results

# save a .RData file with results after computation
save(
  bootstraps,
  file = "data/processed/bootstrapped/bootstrapped-lasso-models.RData"
)
