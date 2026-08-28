# ==============================================================================
# File: R/dgp_misspec_factory.R
# Purpose: Structural-misspecification DGPs for Script 05.
#          Misspecification severity is controlled independently of MIS k.
# ==============================================================================

robust_scale05 <- function(x) {
  s <- stats::mad(x, center = stats::median(x), constant = 1.4826, na.rm = TRUE)
  if (!is.finite(s) || s <= .Machine$double.eps) {
    s <- stats::sd(x, na.rm = TRUE)
  }
  if (!is.finite(s) || s <= .Machine$double.eps) 1 else s
}

# Common severity level -> scenario-specific target.
# Additive misspecification scenarios use robust signal/noise ratios.
# Endogeneity uses target Corr(X, error).
resolve_misspec_severity <- function(
    scenario,
    severity_level,
    additive_grid = c(0, 0.25, 0.50, 1.00, 2.00),
    endogeneity_grid = c(0, 0.10, 0.25, 0.50, 0.75),
    heterosk_grid = c(0, 0.20, 0.50, 1.00, 1.50)
) {
  severity_level <- as.integer(severity_level)
  if (severity_level < 0L || severity_level > 4L) {
    stop("severity_level must be one of 0,1,2,3,4.")
  }
  i <- severity_level + 1L
  
  if (scenario == "correct") return(0)
  if (scenario == "endogeneity") return(endogeneity_grid[i])
  if (scenario == "heteroskedastic") return(heterosk_grid[i])
  additive_grid[i]
}

# Reproduces the distribution families already used in the project, but keeps
# Script 05 self-contained. GPD uses a safe default shape (< 0.5) so mean and
# variance exist in the main experiment; extreme GPD shapes can be supplied
# explicitly as a stress panel.
generate_vector05 <- function(
    n,
    dist_type,
    mix_prop = 0.10,
    skew_t_df = 5,
    gpd_shape = 0.25,
    pareto_shape = 3,
    center = FALSE
) {
  x <- switch(
    dist_type,
    normal = stats::rnorm(n),
    mixed_normal = {
      hi <- stats::runif(n) < mix_prop
      out <- stats::rnorm(n)
      out[hi] <- stats::rnorm(sum(hi), sd = 10)
      out
    },
    skewed_t = {
      if (!requireNamespace("sn", quietly = TRUE)) stop("Install package 'sn'.")
      sn::rst(n, xi = 0, omega = 1, alpha = 5, nu = skew_t_df)
    },
    golm = {
      comp2 <- stats::runif(n) >= 0.7
      out <- stats::rlnorm(n, meanlog = 0, sdlog = 0.5)
      out[comp2] <- stats::rlnorm(sum(comp2), meanlog = 1, sdlog = 1.5)
      out
    },
    beta_logistic = stats::rbeta(n, 2, 5),
    gpd = {
      if (!requireNamespace("evd", quietly = TRUE)) stop("Install package 'evd'.")
      evd::rgpd(n, loc = 0, scale = 1, shape = gpd_shape)
    },
    contaminated = {
      out <- stats::rnorm(n)
      idx <- sample.int(n, size = max(1L, floor(mix_prop * n)))
      out[idx] <- out[idx] + stats::rnorm(length(idx), sd = 50)
      out
    },
    pareto = {
      if (!requireNamespace("actuar", quietly = TRUE)) stop("Install package 'actuar'.")
      actuar::rpareto(n, shape = pareto_shape, scale = 1)
    },
    stop(sprintf("Unsupported distribution '%s'.", dist_type))
  )
  
  if (center) {x <- x - mean(x, na.rm = TRUE)}
  x
}

# Scale an additive misspecification component so its robust scale is a target
# fraction of the error robust scale. Returns the scaled term and implied raw
# coefficient multiplying the unscaled component.
calibrate_additive_misspec <- function(component, error, target_ratio) {
  if (!is.finite(target_ratio) || target_ratio < 0) {
    stop("target_ratio must be finite and non-negative.")
  }
  if (target_ratio == 0) {
    return(list(term = rep(0, length(component)), raw_coef = 0,
                realized_ratio = 0))
  }
  sc_comp <- robust_scale05(component)
  sc_err <- robust_scale05(error)
  raw_coef <- target_ratio * sc_err / sc_comp
  term <- raw_coef * component
  realized <- robust_scale05(term) / sc_err
  list(term = term, raw_coef = raw_coef, realized_ratio = realized)
}

# Generate an auxiliary Z correlated with X while preserving the rank ordering
# of the supplied X distribution. This is used for OVB and interactions.
generate_correlated_z05 <- function(x, rho = 0.5) {
  rho <- max(min(rho, 0.99), -0.99)
  xs <- (x - mean(x)) / stats::sd(x)
  if (any(!is.finite(xs))) xs <- as.numeric(scale(rank(x)))
  z0 <- stats::rnorm(length(x))
  rho * xs + sqrt(1 - rho^2) * z0
}

#' Generate data under controlled structural misspecification.
#'
#' scenario:
#'   correct, ovb, nonlinear, heterogeneous, structural_break,
#'   threshold, missing_interaction, endogeneity, heteroskedastic
#'
#' severity_level = 0,...,4.  It is translated by resolve_misspec_severity().
#' k is deliberately NOT an argument: k belongs to the MIS sensitivity budget,
#' not to general misspecification generation.
generate_misspec_data <- function(
    n,
    x_type = "normal",
    error_type = "normal",
    scenario = "correct",
    severity_level = 0L,
    beta0 = 0,
    beta1 = 1,
    rho_xz = 0.5,
    hetero_group_prop = 0.25,
    break_fraction = 0.75,
    threshold_quantile = 0.75,
    mix_prop = 0.10,
    skew_t_df = 5,
    gpd_shape = 0.25,
    pareto_shape = 3
) {
  allowed <- c(
    "correct", "ovb", "nonlinear", "heterogeneous",
    "structural_break", "threshold", "missing_interaction",
    "endogeneity", "heteroskedastic"
  )
  if (!scenario %in% allowed) stop("Unknown scenario: ", scenario)
  
  x <- generate_vector05(
    n, x_type, mix_prop = mix_prop, skew_t_df = skew_t_df,
    gpd_shape = gpd_shape, pareto_shape = pareto_shape, center = FALSE
  )
  eps0 <- generate_vector05(
    n, error_type, mix_prop = mix_prop, skew_t_df = skew_t_df,
    gpd_shape = gpd_shape, pareto_shape = pareto_shape, center = TRUE
  )
  
  sev <- resolve_misspec_severity(scenario, severity_level)
  z <- generate_correlated_z05(x, rho = rho_xz)
  g <- rep(0L, n)
  d_break <- rep(0L, n)
  miss_component <- rep(0, n)
  miss_term <- rep(0, n)
  raw_coef <- 0
  realized_severity <- 0
  affected_idx <- NULL
  eps <- eps0
  
  wrong_formula <- y ~ x
  correct_formula <- y ~ x
  
  if (scenario == "ovb") {
    miss_component <- z
    cal <- calibrate_additive_misspec(miss_component, eps0, sev)
    miss_term <- cal$term; raw_coef <- cal$raw_coef; realized_severity <- cal$realized_ratio
    correct_formula <- y ~ x + z
    
  } else if (scenario == "nonlinear") {
    # Center X^2 to avoid the misspecification term behaving like an intercept.
    miss_component <- x^2 - stats::median(x^2)
    cal <- calibrate_additive_misspec(miss_component, eps0, sev)
    miss_term <- cal$term; raw_coef <- cal$raw_coef; realized_severity <- cal$realized_ratio
    correct_formula <- y ~ x + I(x^2)
    
  } else if (scenario == "heterogeneous") {
    g <- as.integer(stats::runif(n) < hetero_group_prop)
    miss_component <- x * g
    cal <- calibrate_additive_misspec(miss_component, eps0, sev)
    miss_term <- cal$term; raw_coef <- cal$raw_coef; realized_severity <- cal$realized_ratio
    affected_idx <- if (sev > 0) {
      which(g == 1L)
    } else {
      NULL
    }
    correct_formula <- y ~ x + g + x:g
    
  } else if (scenario == "structural_break") {
    cut <- max(2L, min(n - 1L, floor(n * break_fraction)))
    d_break[(cut + 1L):n] <- 1L
    miss_component <- x * d_break
    cal <- calibrate_additive_misspec(miss_component, eps0, sev)
    miss_term <- cal$term; raw_coef <- cal$raw_coef; realized_severity <- cal$realized_ratio
    affected_idx <- if (sev > 0) {
      which(d_break == 1L)
    } else {
      NULL
    }
    correct_formula <- y ~ x + d_break + x:d_break
    
  } else if (scenario == "threshold") {
    c0 <- as.numeric(stats::quantile(x, probs = threshold_quantile, na.rm = TRUE))
    hinge <- pmax(x - c0, 0)
    miss_component <- hinge
    cal <- calibrate_additive_misspec(miss_component, eps0, sev)
    miss_term <- cal$term; raw_coef <- cal$raw_coef; realized_severity <- cal$realized_ratio
    affected_idx <- which(x > c0)
    correct_formula <- y ~ x + hinge
    
  } else if (scenario == "missing_interaction") {
    miss_component <- x * z
    cal <- calibrate_additive_misspec(miss_component, eps0, sev)
    miss_term <- cal$term; raw_coef <- cal$raw_coef; realized_severity <- cal$realized_ratio
    wrong_formula <- y ~ x + z
    correct_formula <- y ~ x + z + x:z
    
  } else if (scenario == "endogeneity") {
    # Here severity is not an additive signal/noise ratio: it is the target
    # correlation between X and the structural error. Standardize X first.
    xs <- as.numeric(scale(x))
    u <- (eps0 - mean(eps0)) / stats::sd(eps0)
    if (any(!is.finite(u))) u <- stats::rnorm(n)
    eps <- sev * xs + sqrt(max(1 - sev^2, 0)) * u
    eps <- eps * robust_scale05(eps0) / robust_scale05(eps)
    realized_severity <- suppressWarnings(stats::cor(x, eps))
    raw_coef <- sev
    
  } else if (scenario == "heteroskedastic") {
    # Conditional mean remains correct; severity controls the spread gradient.
    xs <- abs(as.numeric(scale(x)))
    mult <- exp(sev * pmin(xs, 3) / 3)
    eps <- eps0 * mult
    eps <- eps - stats::median(eps)
    realized_severity <- sev
    raw_coef <- sev
  }
  
  y <- beta0 + beta1 * x + miss_term + eps
  df <- data.frame(
    y = as.numeric(y), x = as.numeric(x), z = as.numeric(z),
    g = g, d_break = d_break
  )
  if (scenario == "threshold") {
    c0 <- as.numeric(stats::quantile(x, probs = threshold_quantile, na.rm = TRUE))
    df$hinge <- pmax(df$x - c0, 0)
  } else {
    df$hinge <- 0
  }
  
  list(
    data = df,
    scenario = scenario,
    severity_level = severity_level,
    severity_target = sev,
    severity_realized = realized_severity,
    raw_misspec_coef = raw_coef,
    misspec_component = miss_component,
    misspec_term = miss_term,
    error = eps,
    true_beta = beta1,
    wrong_formula = wrong_formula,
    correct_formula = correct_formula,
    affected_idx = affected_idx
  )
}
