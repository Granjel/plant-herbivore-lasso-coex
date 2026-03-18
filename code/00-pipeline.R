# run all scripts

# main analytical pipeline -----------------------------------------------

# set up for analyses (packages, data, etc.)
source("code/01-setup.R")

# empirical interaction coefficients with LV and lasso regularisation
source("code/02-empirical-lasso.R")

# bootstrapped interaction coefficients with LV and lasso regularisation
source("code/03-bootstrapped-lasso.R")
