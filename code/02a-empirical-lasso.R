# empirical interaction coefficients with LV and lasso regularisation

# setup and parallel computation -----------------------------------------

# load general setup and empirical setup
source("code/02-empirical-models.R")

# setup the cluster
n_cores <- parallel::detectCores(logical = TRUE) - 2 # leave 2 cores free
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# assign a distinct random stream to each core
parallel::clusterSetRNGStream(cl, iseed = 1610)

# setup the progress bar
pb <- txtProgressBar(max = length(species), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


# run lasso regularisation models ----------------------------------------

# info and timer start
cat("Running LASSO in parallel on", n_cores, "cores...\n")
start <- Sys.time()

# parallel loop
gli_models <- foreach(
  i = seq_along(species),
  .packages = "glinternet",
  .options.snow = opts
) %dopar%
  {
    # Efficient subsetting
    rows_i <- data$Focal == species[i]

    # The function runs and returns to the list
    glinternet.cv(
      X = data[rows_i, species_subset],
      Y = data[rows_i, "Cover"],
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

# save a .RData file with results after computation, not to run it every time
save(gli_models, file = "data/processed/empirical/empirical-lasso-models.RData")
