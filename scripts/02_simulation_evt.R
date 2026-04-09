# ==============================================================================
# File: /script/02_simulation_evt.R
# Purpose: Execute Phase 2 simulations to evaluate the statistical properties 
#          (size and power) of the EVT-based Most Influential Set hypothesis test.
#          Runs Monte Carlo iterations using the evt_iter() wrapper.
# Inputs: None (Synthetic Data Generation Process defined internally)
# Outputs: ../output/sim_evt_results.rds
# ==============================================================================

library(dplyr)
library(purrr)
source("../R/evt_iter.R")

# 1. Configuration & Parameters
SIM_PARAMS <- list(
  n_iters     = 100,
  n_obs       = 600,
  set_size    = 5,
  block_count = 20,
  seed        = 20260408
)

set.seed(SIM_PARAMS$seed)

# 2. DGP
#' 
#' @param n Number of observations
#' @param set_size Number of observations to randomly select as the "test set"
#' @param inject_influence Logical; if TRUE, forces the selected set to be highly influential
#' @return A list containing y, x, Z, and the indices of the test set
generate_sim_data <- function(n, set_size, inject_influence = FALSE) {
  # Base covariates
  x <- rnorm(n)
  Z <- matrix(rnorm(n * 2), ncol = 2)
  
  # DGP: y = 2x + z1 - 0.5z2 + error
  error <- rnorm(n)
  y <- 2 * x + Z[, 1] - 0.5 * Z[, 2] + error
  
  # Randomly pick a set of indices to test
  test_set <- sample(1:n, size = set_size)
  
  if (isTRUE(inject_influence)) {
    # To engineer influence, we push the test set to extreme leverage (x) 
    # and extreme residual (y) space.
    x[test_set] <- x[test_set] + 5 
    y[test_set] <- y[test_set] - 10 
  }
  
  list(y = y, x = x, Z = Z, set = test_set)
}

# 3. Simulation Loop
cat(sprintf("Starting simulation with %d iterations...\n", SIM_PARAMS$n_iters))

# purrr::map_dfr iterates from 1 to n_iters, runs the block, and row-binds the results
sim_results <- purrr::map_dfr(1:SIM_PARAMS$n_iters, function(i) {
  
  # Progress tracker
  if (i %% 50 == 0) cat(sprintf("  Completed %d / %d...\n", i, SIM_PARAMS$n_iters))
  
  # 3a. Generate Data (Testing the Null: inject_influence = FALSE)
  dat <- generate_sim_data(
    n = SIM_PARAMS$n_obs, 
    set_size = SIM_PARAMS$set_size, 
    inject_influence = FALSE
  )
  
  # Find MIS
  # 1). Fit the OLS model
  df_sim <- data.frame(y = dat$y, x = dat$x, Z1 = dat$Z[,1], Z2 = dat$Z[,2])
  base_model <- lm(y ~ x + Z1 + Z2, data = df_sim)
  
  # 2). Use 'sens' generic.
  sens_obj <- influence::sens(
    base_model,
    lambda = influence::set_lambda("beta_i", pos = 2, sign = sign(coef(base_model)[2]))
  )
  
  # Extract the exact top-k indices of the MIS
  true_mis_indices <- sens_obj$influence$id[1:SIM_PARAMS$set_size]
  
  # 3b. Run the EVT wrapper on the TRUE MIS
  res <- evt_iter(
    y = dat$y, 
    x = dat$x, 
    Z = dat$Z, 
    set = true_mis_indices, 
    block_count = SIM_PARAMS$block_count
  )
  
  # 3c. Append iteration metadata
  res$iter             <- i
  res$n_obs            <- SIM_PARAMS$n_obs
  res$set_size         <- SIM_PARAMS$set_size
  res$inject_influence <- FALSE
  
  return(res)
})

# 4. Save and Quick Diagnostic
if (!dir.exists("../output")) dir.create("../output")

output_file <- "../output/sim_evt_null_distribution.rds"
saveRDS(sim_results, output_file)

cat(sprintf("\nSimulation complete. Saved %d rows to %s\n", nrow(sim_results), output_file))

# Quick sanity check of the convergence rate
convergence_rate <- mean(sim_results$converged, na.rm = TRUE)
cat(sprintf("Optimizer Convergence Rate: %.1f%%\n", convergence_rate * 100))

# Empirical size (False Positive Rate at alpha = 0.05)
# If the math in the paper holds, this should hover right around 0.05
if (convergence_rate > 0) {
  valid_results <- sim_results %>% filter(converged == TRUE)
  empirical_size <- mean(valid_results$p_value < 0.05, na.rm = TRUE)
  cat(sprintf("Empirical Size (alpha = 0.05): %.3f\n", empirical_size))
}
# Under the null, p-values should be uniform on [0,1]
hist(valid_results$p_value, breaks = 20)