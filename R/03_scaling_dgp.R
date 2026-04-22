# ==============================================================================
# File: /R/03_scaling_dgp.R
# Purpose: Data-generating processes for Script 03 (algorithmic comparison).
#          Provides five model architectures of increasing complexity —
#          simple, complex (high-dim + FE), interaction, triple interaction,
#          and nonlinear nuisance (intentional misspecification) — each
#          returning a ready-to-estimate data.frame, formula, and target
#          coefficient position for use with influence::set_lambda().
# Called by: /script/03_algorithmic_comparison.R
# ==============================================================================
#' @param N Integer. Sample size.
#' @param architecture String. One of "simple", "complex", "interaction", 
#'        "triple_interaction", "nonlinear_nuisance".
#' @param beta_target Numeric. The true effect size of the parameter of interest.
#'
#' @return A list containing:
#'   - data: A data.frame with y, X, and relevant Z/M/W controls
#'   - target_pos: Integer indicating the column index of the target parameter in the design matrix
#'   - formula: The formula object to be passed to lm() or the MIS function

generate_scaling_dgp <- function(N, architecture, beta_target = 1.0) {
  
  # Base components present in all models
  X <- rnorm(N, mean = 0, sd = 1)
  epsilon <- rnorm(N, mean = 0, sd = 1)
  
  if (architecture == "simple") {
    Z <- rnorm(N, mean = 0, sd = 1)
    y <- 0.5 + beta_target * X + 1.5 * Z + epsilon
    data <- data.frame(y = y, X = X, Z = Z)
    form <- as.formula(y ~ X + Z)
    target_pos <- 2 # Intercept (1), X (2)
    
  } else if (architecture == "complex") {
    # 20 continuous noise variables
    Z_matrix <- matrix(rnorm(N * 20), nrow = N, ncol = 20)
    # A sparse categorical variable (e.g., 5 states)
    state_fe <- sample(1:5, N, replace = TRUE)
    
    # True DGP only cares about the first 3 Zs and X
    y <- 0.5 + beta_target * X + 1.5 * Z_matrix[,1] - 0.8 * Z_matrix[,2] + epsilon
    
    data <- data.frame(y = y, X = X, Z_matrix, state = as.factor(state_fe))
    form <- as.formula(y ~ X + . - y)
    target_pos <- 2 # Intercept (1), X (2), followed by Zs and FEs
    
  } else if (architecture == "interaction") {
    M <- rnorm(N, mean = 0, sd = 1)
    # y = b0 + b1*X + b2*M + b3*(X*M) + e. Target is the interaction (b3)
    y <- 0.5 + 0.3 * X + 0.4 * M + beta_target * (X * M) + epsilon
    data <- data.frame(y = y, X = X, M = M)
    form <- as.formula(y ~ X * M)
    target_pos <- 4 # Intercept (1), X (2), M (3), X:M (4)
    
  } else if (architecture == "triple_interaction") {
    M <- rnorm(N, mean = 0, sd = 1)
    W <- rnorm(N, mean = 0, sd = 1)
    # Target is the triple interaction (b8)
    # R expands X*M*W to: 1, X, M, W, X:M, X:W, M:W, X:M:W
    y <- 0.5 + 0.2*X + 0.2*M + 0.2*W + 
      0.1*(X*M) + 0.1*(X*W) + 0.1*(M*W) + 
      beta_target * (X * M * W) + epsilon
    data <- data.frame(y = y, X = X, M = M, W = W)
    form <- as.formula(y ~ X * M * W)
    target_pos <- 8 
    
  } else if (architecture == "nonlinear_nuisance") {
    Z <- rnorm(N, mean = 0, sd = 2)
    # True DGP is highly non-linear in Z
    y <- 0.5 + beta_target * X + sin(Z) + 0.5 * (Z^2) + epsilon
    data <- data.frame(y = y, X = X, Z = Z)
    # BUT the estimated model is linear (Misspecification test)
    form <- as.formula(y ~ X + Z)
    target_pos <- 2 
    
  } else {
    stop("Unknown architecture")
  }
  
  return(list(
    data = data,
    target_pos = target_pos,
    formula = form,
    true_beta = beta_target
  ))
}
