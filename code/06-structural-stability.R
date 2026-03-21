# run structural stability analyses on the empirical data

# setup ------------------------------------------------------------------

# load setup
source("code/01-setup.R")

# load empirical matrices
load("data/processed/empirical/empirical-matrices.RData")

# ensure igr is a named vector to match alpha's rownames for subsetting
names(igr) <- rownames(alpha)


# structural coexistence functions ---------------------------------------

# structural niche differences (SND)
Omega <- function(alpha) {
  n <- nrow(alpha)
  Sigma <- solve(t(alpha) %*% alpha, tol = 1e-20)
  d <- pmvnorm(
    lower = rep(0, n),
    upper = rep(Inf, n),
    mean = rep(0, n),
    sigma = Sigma
  )
  out <- log10(d[1]) + n * log10(2)
  return(out)
}

# calculate the centroid of the feasibility domain
r_centroid <- function(alpha) {
  n <- nrow(alpha)
  D <- diag(1 / sqrt(diag(t(alpha) %*% alpha)))
  alpha_n <- alpha %*% D
  r_c <- rowSums(alpha_n) / n
  r_c <- t(t(r_c))
  return(r_c)
}

# structural fitness differences (SFD)
theta <- function(alpha, r) {
  r_c <- r_centroid(alpha)
  out <- acos(sum(r_c * r) / (sqrt(sum(r^2)) * sqrt(sum(r_c^2)))) * 180 / pi
  return(out)
}

# feasibility
test_feasibility <- function(alpha, r) {
  out <- prod(solve(alpha, r) > 0)
  return(out)
}

# intra- and interspecific values, and positive-negative difference
intra <- function(A) sum(diag(A)) / nrow(A)
inter <- function(A) (sum(A) - sum(diag(A))) / (nrow(A) * (nrow(A) - 1))
pnd <- function(A) sum(A)

# compute structural metrics wrapper
structural_coex <- function(alpha, intrinsic, n) {
  combos <- t(combn(rownames(alpha), n))
  results_combos <- as.data.frame(matrix(nrow = dim(combos)[1], ncol = 9))
  row.names(results_combos) <- apply(combos, 1, paste, collapse = "_")
  colnames(results_combos) <- c(
    "SND",
    "SFD",
    "feasibility",
    "intra",
    "inter",
    "iid",
    "skewness",
    "kurtosis",
    "pnd"
  )

  for (i in 1:nrow(combos)) {
    # subset alpha and igr for the specific module
    alpha2 <- as.matrix(alpha[combos[i, ], combos[i, ]])
    intrinsic2 <- as.vector(intrinsic[combos[i, ]])

    # if the interaction matrix is completely empty (no interactions), skip
    if (all(alpha2 == 0)) {
      results_combos[i, ] <- NA
    } else {
      results_combos$SND[i] <- tryCatch(10^Omega(alpha2), error = function(e) {
        NA
      })
      results_combos$SFD[i] <- tryCatch(
        theta(alpha2, intrinsic2),
        error = function(e) NA
      )
      results_combos$feasibility[i] <- tryCatch(
        test_feasibility(alpha2, intrinsic2),
        error = function(e) NA
      )
      results_combos$intra[i] <- tryCatch(intra(alpha2), error = function(e) NA)
      results_combos$inter[i] <- tryCatch(inter(alpha2), error = function(e) NA)
      results_combos$iid[i] <- tryCatch(
        intra(alpha2) - inter(alpha2),
        error = function(e) NA
      )
      results_combos$skewness[i] <- tryCatch(
        skewness(as.numeric(alpha2)),
        error = function(e) NA
      )
      results_combos$kurtosis[i] <- tryCatch(
        kurtosis(as.numeric(alpha2)),
        error = function(e) NA
      )
      results_combos$pnd[i] <- tryCatch(pnd(alpha2), error = function(e) NA)
    }
  }
  results_combos$combos <- rownames(results_combos)
  return(results_combos)
}


# computation ------------------------------------------------------------

# define module size (change between 3 and 4)
richness <- 3

cat("Computing structural metrics for", richness, "-species modules...\n")
res <- structural_coex(alpha = as.matrix(alpha), intrinsic = igr, n = richness)

# clean up results
str_coex <- data.frame(
  "n" = rep(richness, nrow(res)),
  "combos" = res$combos,
  "SND" = res$SND,
  "SFD" = res$SFD,
  "feasibility" = res$feasibility,
  "intra" = res$intra,
  "inter" = res$inter,
  "iid" = res$iid,
  "skewness" = res$skewness,
  "kurtosis" = res$kurtosis,
  "pnd" = res$pnd
)
rownames(str_coex) <- NULL
str_coex <- str_coex[complete.cases(str_coex), ]


# export data and plot ---------------------------------------------------

# save data to results/tables/
write.table(
  str_coex,
  file = paste0("results/tables/str-coex-modules-", richness, "-species.txt"),
  sep = "\t",
  row.names = FALSE
)

# info
cat("Done! Data saved.\n")
