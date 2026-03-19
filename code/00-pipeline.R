# run all scripts

# run computationally heavy analyses?
full_run <- FALSE # TRUE for a full run

# main analytical pipeline -----------------------------------------------

# set up for analyses (packages, data, etc.)
source("code/01-setup.R")

# computationally heavy; run only if full_run is TRUE
if (full_run) {
  # empirical interaction coefficients with LV and lasso regularisation
  source("code/03a-empirical-lasso.R")

  # extract empirical coefficients from glinternet and build matrices
  source("code/03b-empirical-extract-coefs.R")

  # build the matrices for the empirical data
  source("code/03c-empirical-build-matrices.R")

  # bootstrapped interaction coefficients with LV and lasso regularisation
  source("code/04a-bootstrapped-lasso.R")

  # extract bootstrapped coefficients from glinternet and build matrices
  source("code/04b-bootstrapped-extract-coefs.R")

  # build the matrices for the bootstrapped data
  source("code/04c-bootstrapped-build-matrices.R")
}
