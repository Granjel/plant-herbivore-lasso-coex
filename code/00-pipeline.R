# run all scripts

# run computationally heavy analyses?
full_run <- FALSE # TRUE for a full run

# main analytical pipeline -----------------------------------------------

# set up for analyses (packages, data, etc.)
source("code/01-setup.R")

# computationally heavy; run only if full_run is TRUE
if (full_run) {
  # run empirical models
  source("code/02a-empirical-lasso.R")
  source("code/02b-empirical-extract-coefs.R")

  # run bootstrapped models
  source("code/03a-bootstrapped-lasso.R")
  source("code/03b-bootstrapped-extract-coefs.R")
}
