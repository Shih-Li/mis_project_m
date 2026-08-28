# ==============================================================================
# File: R/mis_sensitivity.R
# Purpose: Pure fixed-k coefficient sensitivity for Script 05.
#          Unlike run_mis_directional() in Script 04, direction is selected by
#          maximum absolute movement away from the full OLS coefficient.
# ==============================================================================

fit_target_ols05 <- function(formula, data, target_var = "x", exclude_idx = integer(0)) {
  keep <- if (length(exclude_idx)) setdiff(seq_len(nrow(data)), exclude_idx) else seq_len(nrow(data))
  tryCatch({
    mod <- stats::lm(formula, data = data[keep, , drop = FALSE])
    sm <- summary(mod)$coefficients
    if (!target_var %in% rownames(sm)) return(c(coef = NA_real_, se = NA_real_))
    c(coef = unname(sm[target_var, "Estimate"]), se = unname(sm[target_var, "Std. Error"]))
  }, error = function(e) c(coef = NA_real_, se = NA_real_))
}

fit_target_mm05 <- function(formula, data, target_var = "x", exclude_idx = integer(0)) {
  if (!requireNamespace("robustbase", quietly = TRUE)) {
    return(c(coef = NA_real_, se = NA_real_))
  }
  keep <- if (length(exclude_idx)) setdiff(seq_len(nrow(data)), exclude_idx) else seq_len(nrow(data))
  tryCatch({
    mod <- robustbase::lmrob(formula, data = data[keep, , drop = FALSE], setting = "KS2014")
    sm <- summary(mod)$coefficients
    if (!target_var %in% rownames(sm)) return(c(coef = NA_real_, se = NA_real_))
    c(coef = unname(sm[target_var, "Estimate"]), se = unname(sm[target_var, "Std. Error"]))
  }, error = function(e) c(coef = NA_real_, se = NA_real_))
}

# Optional Rcpp loader. If compilation fails, all functions transparently fall
# back to the existing pure-R dinkelbach_topk_lm().
load_mis_cpp_kernel <- function(cpp_file = "../src/dinkelbach_topk_cpp.cpp", quiet = TRUE) {
  if (!requireNamespace("Rcpp", quietly = TRUE) || !file.exists(cpp_file)) return(FALSE)
  ok <- tryCatch({
    Rcpp::sourceCpp(cpp_file, rebuild = FALSE, showOutput = !quiet, verbose = FALSE)
    exists("dinkelbach_topk_cpp", mode = "function")
  }, error = function(e) FALSE)
  isTRUE(ok)
}

# FWL wrapper that uses the optional C++ low-level Dinkelbach kernel.
dinkelbach_topk_lm05 <- function(mod, pos = 2L, sign = 1L, k = 1L, use_cpp = TRUE) {
  if (!use_cpp || !exists("dinkelbach_topk_cpp", mode = "function")) {
    return(dinkelbach_topk_lm(mod, pos = pos, sign = sign, k = k))
  }
  
  X <- stats::model.matrix(mod)
  y <- stats::model.response(stats::model.frame(mod))
  p <- ncol(X)
  if (pos < 1L || pos > p) stop("Invalid target position.")
  
  if (p == 1L) {
    x_fwl <- X[, 1L]
    r_fwl <- stats::residuals(mod)
  } else {
    Z <- X[, setdiff(seq_len(p), pos), drop = FALSE]
    qz <- qr(Z)
    x_fwl <- qr.resid(qz, X[, pos])
    y_fwl <- qr.resid(qz, y)
    bj <- sum(x_fwl * y_fwl) / sum(x_fwl^2)
    r_fwl <- y_fwl - x_fwl * bj
  }
  
  ans <- dinkelbach_topk_cpp(
    x = as.numeric(x_fwl), r = as.numeric(r_fwl), k = as.integer(k),
    sgn = as.integer(sign), sum_x2 = sum(x_fwl^2), max_iter = 50L, tol = 1e-9
  )
  as.integer(ans$indices)
}

# ==============================================================================
# Prepare MIS problem once per fitted OLS model
# ==============================================================================

prepare_mis_problem05 <- function(
    mod,
    pos = 2L
) {
  
  X <- stats::model.matrix(mod)
  
  y <- stats::model.response(
    stats::model.frame(mod)
  )
  
  p <- ncol(X)
  
  
  if (pos < 1L || pos > p) {
    stop("Invalid target coefficient position.")
  }
  
  
  # --------------------------------------------------------------
  # FWL residualisation
  # --------------------------------------------------------------
  
  if (p == 1L) {
    
    x_fwl <- X[, 1L]
    
    r_fwl <- stats::residuals(mod)
    
  } else {
    
    Z_cols <- setdiff(
      seq_len(p),
      pos
    )
    
    
    Z <- X[
      ,
      Z_cols,
      drop = FALSE
    ]
    
    
    qz <- qr(Z)
    
    
    x_fwl <- qr.resid(
      qz,
      X[, pos]
    )
    
    
    y_fwl <- qr.resid(
      qz,
      y
    )
    
    
    beta_fwl <-
      sum(x_fwl * y_fwl) /
      sum(x_fwl^2)
    
    
    r_fwl <-
      y_fwl -
      x_fwl * beta_fwl
  }
  
  
  list(
    
    x_fwl =
      as.numeric(x_fwl),
    
    r_fwl =
      as.numeric(r_fwl),
    
    sum_x2 =
      sum(x_fwl^2)
    
  )
}


# ==============================================================================
# Run Dinkelbach from an already prepared FWL problem
# ==============================================================================

run_dinkelbach_problem05 <- function(
    problem,
    sign,
    k,
    use_cpp = TRUE
) {
  
  k <- as.integer(k)
  
  
  # --------------------------------------------------------------
  # C++ kernel
  # --------------------------------------------------------------
  
  if (
    isTRUE(use_cpp) &&
    exists(
      "dinkelbach_topk_cpp",
      mode = "function"
    )
  ) {
    
    ans <- dinkelbach_topk_cpp(
      
      x =
        problem$x_fwl,
      
      r =
        problem$r_fwl,
      
      k =
        k,
      
      sgn =
        as.integer(sign),
      
      sum_x2 =
        problem$sum_x2,
      
      max_iter =
        50L,
      
      tol =
        1e-9
    )
    
    
    return(
      as.integer(
        ans$indices
      )
    )
  }
  
  
  # --------------------------------------------------------------
  # Pure-R fallback
  # --------------------------------------------------------------
  
  ans <- dinkelbach_topk(
    
    x =
      problem$x_fwl,
    
    r =
      problem$r_fwl,
    
    k =
      k,
    
    sgn =
      as.integer(sign),
    
    sum_x2 =
      problem$sum_x2
    
  )
  
  
  as.integer(
    ans$indices
  )
}

run_mis_sensitivity <- function(
    
  mod_full,
  
  formula,
  
  data,
  
  k,
  
  target_var = "x",
  
  target_pos = 2L,
  
  use_cpp = TRUE,
  
  full_fit = NULL,
  
  problem = NULL
  
) {
  
  k <- as.integer(k)
  
  
  # --------------------------------------------------------------
  # Reuse full OLS if already computed
  # --------------------------------------------------------------
  
  if (is.null(full_fit)) {
    
    full_fit <- fit_target_ols05(
      formula,
      data,
      target_var
    )
  }
  
  
  if (
    k <= 0L ||
    !is.finite(
      full_fit["coef"]
    )
  ) {
    
    return(
      list(
        
        indices =
          integer(0),
        
        direction =
          0L,
        
        full =
          full_fit,
        
        deleted =
          full_fit,
        
        delta_beta =
          0,
        
        standardized_delta =
          0
      )
    )
  }
  
  
  # --------------------------------------------------------------
  # Reuse FWL problem if already prepared
  # --------------------------------------------------------------
  
  if (is.null(problem)) {
    
    problem <- prepare_mis_problem05(
      mod_full,
      pos = target_pos
    )
  }
  
  
  # --------------------------------------------------------------
  # Positive MIS direction
  # --------------------------------------------------------------
  
  idx_pos <- run_dinkelbach_problem05(
    
    problem =
      problem,
    
    sign =
      1L,
    
    k =
      k,
    
    use_cpp =
      use_cpp
  )
  
  
  # --------------------------------------------------------------
  # Negative MIS direction
  # --------------------------------------------------------------
  
  idx_neg <- run_dinkelbach_problem05(
    
    problem =
      problem,
    
    sign =
      -1L,
    
    k =
      k,
    
    use_cpp =
      use_cpp
  )
  
  
  # --------------------------------------------------------------
  # Only two OLS deletion fits are required
  # --------------------------------------------------------------
  
  pos <- fit_target_ols05(
    
    formula,
    data,
    target_var,
    idx_pos
  )
  
  
  neg <- fit_target_ols05(
    
    formula,
    data,
    target_var,
    idx_neg
  )
  
  
  dpos <- abs(
    pos["coef"] -
      full_fit["coef"]
  )
  
  
  dneg <- abs(
    neg["coef"] -
      full_fit["coef"]
  )
  
  
  choose_pos <-
    is.finite(dpos) &&
    (
      !is.finite(dneg) ||
        dpos >= dneg
    )
  
  
  if (choose_pos) {
    
    idx <- idx_pos
    
    deleted <- pos
    
    direction <- 1L
    
  } else {
    
    idx <- idx_neg
    
    deleted <- neg
    
    direction <- -1L
  }
  
  
  delta <- abs(
    deleted["coef"] -
      full_fit["coef"]
  )
  
  
  standardized_delta <-
    
    if (
      is.finite(full_fit["se"]) &&
      full_fit["se"] > 0
    ) {
      
      delta /
        full_fit["se"]
      
    } else {
      
      NA_real_
      
    }
  
  
  list(
    
    indices =
      idx,
    
    direction =
      direction,
    
    full =
      full_fit,
    
    deleted =
      deleted,
    
    delta_beta =
      unname(delta),
    
    standardized_delta =
      unname(
        standardized_delta
      )
  )
}

refit_after_selection05 <- function(data, formula, indices, target_var = "x") {
  list(
    ols = fit_target_ols05(formula, data, target_var, indices),
    mm  = fit_target_mm05(formula, data, target_var, indices)
  )
}
