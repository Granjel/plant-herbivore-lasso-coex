# bootstrapped interaction coefficients with LV and lasso regularisation

# setup and compute bootstrapped data ------------------------------------

# load lasso parameters
source("code/02-lasso-parameters.R")

# times to resample (number of bootstraps)
t_boot <- 500

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
n_cores <- parallel::detectCores(logical = TRUE) - 2 # leave 2 cores free
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# assign a distinct random stream to each core for perfect reproducibility
parallel::clusterSetRNGStream(cl, iseed = 1610)

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
    # Extract the specific resampled dataframe for this iteration
    dat_i <- bootstraps$resampled_data[[i]]

    # The function runs and returns the model to the list
    glinternet.cv(
      X = dat_i[, species_subset],
      Y = dat_i$Cover,
      numLevels = rep(1, n_col),
      nLambda = max_lambda
    )
  }

# stop the cluster and close the progress bar
close(pb)
stopCluster(cl)

# timer end
end <- Sys.time()
print(end - start)

# sound notification when done
beep(1)


# save results -----------------------------------------------------------

# bind the results back to the tibble
bootstraps$gli_models <- lasso_results

# save a .RData file with results after computation
save(
  bootstraps,
  file = "data/processed/bootstrapped/bootstrapped-lasso-models.RData"
)
