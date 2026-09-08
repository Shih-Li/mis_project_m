# ==============================================================================
# File: /scripts/05a_detector.R
#
# FORMAL 05a -- MIS AS A MISSPECIFICATION DETECTOR
#
# This is a NEW experiment. It does not modify, re-run, or overwrite anything
# belonging to script 05. Frozen 05 remains the estimation/robustness evidence
# (bias, RMSE, coverage, MM comparisons). 05a estimates the operating
# characteristic of MIS as a detector: size, power, and detection boundary.
#
# ------------------------------------------------------------------------------
# DESIGN DECISIONS, AND THE PILOT EVIDENCE BEHIND THEM
# ------------------------------------------------------------------------------
#
# 1. SEVERITY IS SD-CALIBRATED:  lambda = sd(term) / sd(eps0)
#
#    Frozen 05 used robust_scale05 (MAD). MAD only sees the bulk of the
#    component, so for x^2 under contaminated X the same nominal severity
#    delivered a physical signal ~1000x larger than under normal X. Pilot 3
#    measured a 620x boundary ratio across X types under MAD versus 2.81x
#    under SD; pilot 4 confirmed SD over a quantile-based alternative
#    (total boundary spread 3.74x versus 6.68x) and showed SD is invariant
#    across error distributions (boundary 0.241 with normal errors, 0.247
#    with golm errors, at normal X).
#
#    NOTE: robust_scale05() silently falls back to sd() when MAD is ~0. The
#    threshold hinge is zero for ~75% of observations, so frozen 05 was
#    ALREADY sd-calibrating that one scenario. SD here makes the treatment
#    uniform rather than accidental.
#
# 2. CALIBRATION IS INDEPENDENT OF EVALUATION
#
#    Frozen 05 took the 95% cutoff from the same severity-0 draws it then
#    scored, under shared CRN seeds -- so its null rejection rate of exactly
#    0.05 was arithmetic, not evidence. Here a disjoint seed stream produces
#    the cutoffs. CRN across severity is retained WITHIN the evaluation
#    sample, which is where it belongs.
#
# 3. CALIBRATION GROUPS ON NULL MODEL CLASS, NOT SCENARIO
#
#    At severity 0 every scenario collapses to the same DGP (the endogeneity
#    branch reduces algebraically to eps0 exactly). But missing_interaction
#    fits y ~ x + z, so its FWL residualisation and null distribution differ.
#    Two classes, not nine: "x" and "x_plus_z".
#
# 4. THE NULL HAS NO FIXED DISTRIBUTION
#
#    Pilot 3: log-log slope of null MIS on n is 0.493 (normal) and 0.480
#    (contaminated), R^2 > 0.999. Standardized delta diverges as sqrt(n).
#    Per-environment empirical calibration is therefore mandatory.
#
# 5. COMPARATORS ARE SIZE-CORRECTED
#
#    Pilot 1-2: RESET ran at nominal size 0.085-0.110 even with vcovHC.
#    Comparing a size-0.10 test with a size-0.05 MIS is not a contest. The
#    calibration stream supplies the empirical 5th percentile of each test's
#    null p-value distribution as its critical value.
#
# 6. NO lmrob ANYWHERE
#
#    Frozen 05 spent ~98% of its 108 wall-hours inside lmrob. None of it is
#    detector evidence. 05a is OLS-only.
#
# 7. structural_break IS DROPPED
#
#    With i.i.d. rows, "last 25% of row numbers get another slope" is
#    distributionally identical to the heterogeneous subgroup scenario, which
#    is why frozen 05 reported 17.868 vs 17.881 and power 0.529 vs 0.537.
#    Its budget is reallocated to endogeneity_nl, which tests whether the
#    endogeneity invariance result is specific to linear-in-x confounding.
#
# ------------------------------------------------------------------------------
# COST
# ------------------------------------------------------------------------------
#   Calibration : 192 cells x 2500 draws   ~ 2.5 h wall on 12 workers
#   Evaluation  : 7680 cells x  500 draws  ~ 20  h wall on 12 workers
#   Set n_iters_eval = 1000 to halve the Monte Carlo error at double the time.
#
# ------------------------------------------------------------------------------
# OUTPUT
# ------------------------------------------------------------------------------
#   ../output/temp/05a_detector/cal_cell_XXXXXX.rds     (192 files)
#   ../output/temp/05a_detector/eval_cell_XXXXXX.rds    (7680 files)
#   ../output/05a_null_cutoffs.rds
#   ../output/05a_design_grids.rds
#
# Post-processing (size, power, isotonic boundary) lives in a separate script
# so that this one can be re-entered safely at any point.
# ==============================================================================


rm(list = ls())


# ==============================================================================
# 1. Libraries
# ==============================================================================

suppressPackageStartupMessages({
  library(future)
  library(furrr)
  library(parallel)
})

stopifnot(requireNamespace("lmtest", quietly = TRUE))
stopifnot(requireNamespace("sandwich", quietly = TRUE))


# ==============================================================================
# 2. Frozen modules -- READ ONLY. Nothing below modifies these files.
# ==============================================================================

source("../R/helpers_local.R")
source("../R/utils_checkpoint.R")
source("../R/diagnostics_classical.R")
source("../R/dinkelbach_topk.R")
source("../R/dgp_misspec_factory.R")
source("../R/mis_sensitivity.R")


# ==============================================================================
# 3. Parameters
# ==============================================================================

sim_params05a <- list(
  
  n_iters_cal  = 2500L,
  n_iters_eval = 500L,
  
  # Disjoint seed streams. Evaluation shares 05's base so that environments
  # remain comparable; calibration is offset far enough to guarantee no
  # collision with any 05 or 05a evaluation seed.
  seed_eval = 20260828L,
  seed_cal  = 77260828L,
  
  use_cpp = TRUE,
  
  workers = max(1L, min(future::availableCores() - 1L, 12L)),
  
  gpd_shape    = 0.25,
  pareto_shape = 3,
  skew_t_df    = 5,
  
  mix_prop = 0.10,
  
  rho_xz            = 0.50,
  hetero_group_prop = 0.25,
  threshold_quantile = 0.75,
  
  alpha = 0.05
)


# ==============================================================================
# 4. Grids
# ==============================================================================

n_grid05a <- c(500L, 1000L, 2500L, 5000L)

x_grid05a <- c("normal", "mixed_normal", "contaminated")

error_grid05a <- c("normal", "mixed_normal", "skewed_t", "golm",
                   "beta_logistic", "gpd", "contaminated", "pareto")

k_fraction_grid05a <- c(0.010, 0.025, 0.050, 0.100)

scenario_grid05a <- c("nonlinear", "threshold", "heteroskedastic",
                      "missing_interaction", "ovb", "heterogeneous",
                      "endogeneity", "endogeneity_nl")

# ------------------------------------------------------------------------------
# Severity grids. Ten points each, in the units named.
#
# additive / threshold : lambda = sd(term)/sd(eps0)
# heteroskedastic      : spread-gradient parameter (multiplicative, unchanged)
# endogeneity          : target Corr(X, error)
#
# The additive top of 1.40 is set by the hardest cell pilot 4 found:
# contaminated X with golm errors went 0.183 -> 1.000 in power between
# lambda 0.60 and 1.00, boundary 0.902.
# ------------------------------------------------------------------------------

severity_grids05a <- list(
  additive        = c(0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.45, 0.65, 0.95, 1.40),
  threshold       = c(0, 0.04, 0.08, 0.12, 0.16, 0.24, 0.32, 0.48, 0.64, 1.00),
  heteroskedastic = c(0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.55, 0.75, 1.00),
  endogeneity     = c(0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.45, 0.60, 0.75, 0.90)
)

severity_family05a <- function(scenario) {
  if (scenario == "threshold") return("threshold")
  if (scenario == "heteroskedastic") return("heteroskedastic")
  if (scenario %in% c("endogeneity", "endogeneity_nl")) return("endogeneity")
  "additive"
}

# Correct-state (respecification) rows are expensive and are not needed at
# every severity. Indices into the 10-point grid.
correct_state_idx05a <- c(1L, 4L, 7L, 10L)

# Scenarios with a meaningful correct specification.
has_correct05a <- c("nonlinear", "threshold", "missing_interaction",
                    "ovb", "heterogeneous")

# Null model class implied by each scenario's WRONG formula.
null_class05a <- function(scenario) {
  if (scenario == "missing_interaction") "x_plus_z" else "x"
}

null_class_grid05a <- c("x", "x_plus_z")


# ==============================================================================
# 5. 05a-local functions
#
# Deliberately defined here rather than by editing ../R/*. Frozen 05 shares
# those files; a "backward-compatible" edit is still a mutation. These reuse
# the frozen primitives (generate_vector05, generate_correlated_z05,
# robust_scale05) but implement the 05a severity definition.
# ==============================================================================

# ------------------------------------------------------------------------------
# SD-based additive severity calibration (replaces calibrate_additive_misspec)
# ------------------------------------------------------------------------------

calibrate_sd05a <- function(component, error, target_ratio) {
  
  if (!is.finite(target_ratio) || target_ratio <= 0) {
    return(list(term = rep(0, length(component)), raw_coef = 0,
                realized_ratio = 0))
  }
  
  sc_comp <- stats::sd(component, na.rm = TRUE)
  sc_err  <- stats::sd(error, na.rm = TRUE)
  
  if (!is.finite(sc_comp) || sc_comp <= .Machine$double.eps) {
    return(list(term = rep(0, length(component)), raw_coef = 0,
                realized_ratio = 0))
  }
  
  raw_coef <- target_ratio * sc_err / sc_comp
  term     <- raw_coef * component
  
  list(term = term, raw_coef = raw_coef,
       realized_ratio = stats::sd(term, na.rm = TRUE) / sc_err)
}


# ------------------------------------------------------------------------------
# Population 75th percentile of X, per x_type.
#
# qnorm(.75) is correct ONLY for x_type = "normal". Using it everywhere would
# let the affected fraction drift with the X distribution.
# ------------------------------------------------------------------------------

pop_q75_05a <- function(x_type, mix_prop = 0.10) {
  if (x_type == "normal") return(stats::qnorm(0.75))
  if (x_type == "mixed_normal")
    return(stats::uniroot(function(q) (1 - mix_prop) * stats::pnorm(q) +
                            mix_prop * stats::pnorm(q / 10) - 0.75,
                          c(-20, 20), tol = 1e-10)$root)
  if (x_type == "contaminated")
    return(stats::uniroot(function(q) (1 - mix_prop) * stats::pnorm(q) +
                            mix_prop * stats::pnorm(q / sqrt(1 + 50^2)) - 0.75,
                          c(-10, 10), tol = 1e-10)$root)
  stop("pop_q75_05a: unsupported x_type '", x_type, "'")
}

c0_pop05a <- vapply(x_grid05a, pop_q75_05a, numeric(1),
                    mix_prop = sim_params05a$mix_prop)


# ------------------------------------------------------------------------------
# DGP
#
# RNG consumption is held constant across severity within a scenario, so the
# CRN seed reproduces identical X, eps0, z and g while only the severity
# parameter moves.
# ------------------------------------------------------------------------------

generate_misspec_data05a <- function(n, x_type, error_type, scenario, sev,
                                     c0_pop, mix_prop = 0.10, rho_xz = 0.50,
                                     hetero_group_prop = 0.25,
                                     beta0 = 0, beta1 = 1,
                                     gpd_shape = 0.25, pareto_shape = 3,
                                     skew_t_df = 5) {
  
  x <- generate_vector05(n, x_type, mix_prop = mix_prop,
                         skew_t_df = skew_t_df, gpd_shape = gpd_shape,
                         pareto_shape = pareto_shape, center = FALSE)
  
  eps0 <- generate_vector05(n, error_type, mix_prop = mix_prop,
                            skew_t_df = skew_t_df, gpd_shape = gpd_shape,
                            pareto_shape = pareto_shape, center = TRUE)
  
  z <- generate_correlated_z05(x, rho = rho_xz)
  g <- as.integer(stats::runif(n) < hetero_group_prop)
  
  miss_term <- rep(0, n)
  eps       <- eps0
  raw_coef  <- 0
  realized  <- 0
  affected  <- NULL
  hinge     <- rep(0, n)
  
  wrong_formula   <- y ~ x
  correct_formula <- y ~ x
  
  if (scenario == "nonlinear") {
    cal <- calibrate_sd05a(x^2 - stats::median(x^2), eps0, sev)
    miss_term <- cal$term; raw_coef <- cal$raw_coef; realized <- cal$realized_ratio
    correct_formula <- y ~ x + I(x^2)
    
  } else if (scenario == "threshold") {
    hinge <- pmax(x - c0_pop, 0)
    cal <- calibrate_sd05a(hinge, eps0, sev)
    miss_term <- cal$term; raw_coef <- cal$raw_coef; realized <- cal$realized_ratio
    if (sev > 0) affected <- which(x > c0_pop)
    correct_formula <- y ~ x + hinge
    
  } else if (scenario == "ovb") {
    cal <- calibrate_sd05a(z, eps0, sev)
    miss_term <- cal$term; raw_coef <- cal$raw_coef; realized <- cal$realized_ratio
    correct_formula <- y ~ x + z
    
  } else if (scenario == "heterogeneous") {
    cal <- calibrate_sd05a(x * g, eps0, sev)
    miss_term <- cal$term; raw_coef <- cal$raw_coef; realized <- cal$realized_ratio
    if (sev > 0) affected <- which(g == 1L)
    correct_formula <- y ~ x + g + x:g
    
  } else if (scenario == "missing_interaction") {
    cal <- calibrate_sd05a(x * z, eps0, sev)
    miss_term <- cal$term; raw_coef <- cal$raw_coef; realized <- cal$realized_ratio
    wrong_formula   <- y ~ x + z
    correct_formula <- y ~ x + z + x:z
    
  } else if (scenario == "endogeneity") {
    xs <- as.numeric(scale(x))
    u  <- (eps0 - mean(eps0)) / stats::sd(eps0)
    if (any(!is.finite(u))) u <- stats::rnorm(n)
    eps <- sev * xs + sqrt(max(1 - sev^2, 0)) * u
    eps <- eps * robust_scale05(eps0) / robust_scale05(eps)
    realized <- suppressWarnings(stats::cor(x, eps))
    raw_coef <- sev
    
  } else if (scenario == "endogeneity_nl") {
    # Nonlinear confounding. The linear-endogeneity invariance proof relies on
    # the confounder being a linear function of x; squaring breaks that, so MIS
    # should become detectable here. This is the scope test for that result.
    xs <- as.numeric(scale(x))
    q  <- xs^2 - mean(xs^2)
    qs <- if (stats::sd(q) > 0) q / stats::sd(q) else q
    u  <- (eps0 - mean(eps0)) / stats::sd(eps0)
    if (any(!is.finite(u))) u <- stats::rnorm(n)
    eps <- sev * qs + sqrt(max(1 - sev^2, 0)) * u
    eps <- eps * robust_scale05(eps0) / robust_scale05(eps)
    realized <- suppressWarnings(stats::cor(x^2, eps))
    raw_coef <- sev
    
  } else if (scenario == "heteroskedastic") {
    xs   <- abs(as.numeric(scale(x)))
    eps  <- eps0 * exp(sev * pmin(xs, 3) / 3)
    eps  <- eps - mean(eps, na.rm = TRUE)
    realized <- sev
    raw_coef <- sev
    
  } else if (scenario != "correct") {
    stop("generate_misspec_data05a: unknown scenario '", scenario, "'")
  }
  
  df <- data.frame(y = as.numeric(beta0 + beta1 * x + miss_term + eps),
                   x = as.numeric(x), z = as.numeric(z),
                   g = g, hinge = as.numeric(hinge))
  
  list(data = df, severity_target = sev, severity_realized = realized,
       raw_misspec_coef = raw_coef, wrong_formula = wrong_formula,
       correct_formula = correct_formula, affected_idx = affected)
}


# ------------------------------------------------------------------------------
# Specification tests. Auxiliary regressions are built from the model's own
# regressors so the same code serves both null model classes.
#
# BP uses ~ x + I(x^2) rather than ~ x: pilot 2 showed the heteroskedastic DGP
# is U-shaped in |x|, where a linear auxiliary regression is weak
# (power 0.505 vs 1.000 at maximum severity, n = 5000).
# ------------------------------------------------------------------------------

spec_tests05a <- function(mod, dat, regressors) {
  
  vf <- stats::as.formula(
    paste("~", paste(c(regressors, sprintf("I(%s^2)", regressors)),
                     collapse = " + ")))
  
  reset_p <- tryCatch(
    lmtest::resettest(mod, power = 2:3, type = "regressor",
                      vcov = sandwich::vcovHC)$p.value,
    error = function(e) NA_real_)
  
  bp_p <- tryCatch(
    lmtest::bptest(mod, varformula = vf, studentize = TRUE,
                   data = dat)$p.value,
    error = function(e) NA_real_)
  
  c(reset_p = unname(reset_p), bp_p = unname(bp_p))
}


# ------------------------------------------------------------------------------
# All influence statistics for one fitted model, at every k, plus the two
# specification tests (which do not depend on k).
# ------------------------------------------------------------------------------

fit_stats05a <- function(formula, dat, target_var, k_fractions, use_cpp,
                         affected_idx = NULL, do_spec_tests = TRUE) {
  
  n <- nrow(dat)
  
  mod <- tryCatch(stats::lm(formula, data = dat), error = function(e) NULL)
  if (is.null(mod)) return(NULL)
  
  sm <- tryCatch(summary(mod)$coefficients, error = function(e) NULL)
  if (is.null(sm) || !(target_var %in% rownames(sm))) return(NULL)
  
  b <- unname(sm[target_var, "Estimate"])
  s <- unname(sm[target_var, "Std. Error"])
  if (!is.finite(b) || !is.finite(s) || s <= 0) return(NULL)
  
  target_pos <- match(target_var, rownames(sm))
  
  cl <- tryCatch(get_all_classical(mod, target_var), error = function(e) NULL)
  
  rk <- if (is.null(cl)) NULL else list(
    Cook     = cl$id[order(cl$cooks_d, decreasing = TRUE)],
    DFBETAS  = cl$id[order(abs(cl$dfbetas_target), decreasing = TRUE)],
    Leverage = cl$id[order(cl$leverage, decreasing = TRUE)])
  
  rows <- vector("list", length(k_fractions))
  
  for (j in seq_along(k_fractions)) {
    
    kf <- k_fractions[j]
    k  <- max(1L, as.integer(round(kf * n)))
    
    m <- tryCatch(
      run_mis_sensitivity(mod_full = mod, formula = formula, data = dat,
                          k = k, target_var = target_var,
                          target_pos = target_pos, use_cpp = use_cpp),
      error = function(e) NULL)
    
    mis_stat <- if (is.null(m)) NA_real_ else m$standardized_delta
    mis_idx  <- if (is.null(m)) integer(0) else m$indices
    
    cls <- c(Cook = NA_real_, DFBETAS = NA_real_, Leverage = NA_real_)
    if (!is.null(rk)) {
      for (nm in names(rk)) {
        idx <- rk[[nm]][seq_len(k)]
        f <- tryCatch(fit_target_ols05(formula, dat, target_var, idx),
                      error = function(e) NULL)
        if (!is.null(f)) cls[nm] <- unname(abs(f["coef"] - b) / s)
      }
    }
    
    # Subgroup recovery for MIS only. Recall is bounded above by
    # k / |affected|, so precision and lift are the interpretable columns.
    prec <- NA_real_; rec <- NA_real_; lift <- NA_real_; afrac <- NA_real_
    if (!is.null(affected_idx) && length(affected_idx) > 0 &&
        length(mis_idx) > 0) {
      hit   <- sum(mis_idx %in% affected_idx)
      prec  <- hit / length(mis_idx)
      rec   <- hit / length(affected_idx)
      afrac <- length(affected_idx) / n
      lift  <- if (afrac > 0) prec / afrac else NA_real_
    }
    
    rows[[j]] <- data.frame(
      k_fraction = kf, k = k,
      full_coef = b, full_se = s,
      MIS = mis_stat, Cook = unname(cls["Cook"]),
      DFBETAS = unname(cls["DFBETAS"]), Leverage = unname(cls["Leverage"]),
      mis_precision = prec, mis_recall = rec, mis_lift = lift,
      affected_fraction = afrac,
      stringsAsFactors = FALSE)
  }
  
  out <- do.call(rbind, rows)
  
  if (do_spec_tests) {
    regs <- attr(stats::terms(mod), "term.labels")
    regs <- regs[!grepl(":", regs, fixed = TRUE)]
    st <- spec_tests05a(mod, dat, regs)
    out$reset_p <- unname(st["reset_p"])
    out$bp_p    <- unname(st["bp_p"])
  } else {
    out$reset_p <- NA_real_
    out$bp_p    <- NA_real_
  }
  
  out
}


# ------------------------------------------------------------------------------
# One CALIBRATION draw. Severity 0, null model class fixed.
# ------------------------------------------------------------------------------

run_cal_iter05a <- function(iter, n, x_type, error_type, null_class,
                            k_fractions, seed_base, use_cpp, params, c0) {
  
  set.seed(seed_base + iter)
  
  d <- tryCatch(
    generate_misspec_data05a(
      n = n, x_type = x_type, error_type = error_type,
      scenario = "correct", sev = 0, c0_pop = c0,
      mix_prop = params$mix_prop, rho_xz = params$rho_xz,
      hetero_group_prop = params$hetero_group_prop,
      gpd_shape = params$gpd_shape, pareto_shape = params$pareto_shape,
      skew_t_df = params$skew_t_df),
    error = function(e) NULL)
  
  if (is.null(d)) return(NULL)
  
  fml <- if (null_class == "x_plus_z") y ~ x + z else y ~ x
  
  st <- fit_stats05a(fml, d$data, "x", k_fractions, use_cpp,
                     affected_idx = NULL, do_spec_tests = TRUE)
  
  if (is.null(st)) return(NULL)
  
  st$iter <- iter
  st
}


# ------------------------------------------------------------------------------
# One EVALUATION draw. Wrong state always; correct state only on the reduced
# severity subgrid, and only for scenarios that have a correct specification.
# ------------------------------------------------------------------------------

run_eval_iter05a <- function(iter, n, x_type, error_type, scenario,
                             severity_index, severity_target, k_fractions,
                             seed_base, use_cpp, params, c0,
                             run_correct) {
  
  set.seed(seed_base + iter)
  
  d <- tryCatch(
    generate_misspec_data05a(
      n = n, x_type = x_type, error_type = error_type,
      scenario = scenario, sev = severity_target, c0_pop = c0,
      mix_prop = params$mix_prop, rho_xz = params$rho_xz,
      hetero_group_prop = params$hetero_group_prop,
      gpd_shape = params$gpd_shape, pareto_shape = params$pareto_shape,
      skew_t_df = params$skew_t_df),
    error = function(e) NULL)
  
  if (is.null(d)) return(NULL)
  
  out <- list()
  
  wrong <- fit_stats05a(d$wrong_formula, d$data, "x", k_fractions, use_cpp,
                        affected_idx = d$affected_idx, do_spec_tests = TRUE)
  
  if (!is.null(wrong)) {
    wrong$model_state <- "wrong"
    out[[length(out) + 1L]] <- wrong
  }
  
  if (isTRUE(run_correct)) {
    corr <- fit_stats05a(d$correct_formula, d$data, "x", k_fractions, use_cpp,
                         affected_idx = d$affected_idx, do_spec_tests = FALSE)
    if (!is.null(corr)) {
      corr$model_state <- "correct"
      out[[length(out) + 1L]] <- corr
    }
  }
  
  if (length(out) == 0L) return(NULL)
  
  res <- do.call(rbind, out)
  res$iter <- iter
  res$severity_realized <- d$severity_realized
  res$raw_misspec_coef  <- d$raw_misspec_coef
  res
}


# ==============================================================================
# 6. Design grids
# ==============================================================================

# ------------------------------------------------------------------------------
# Calibration: 96 environments x 2 null model classes. No scenario axis and no
# severity axis, because severity 0 is one DGP.
# ------------------------------------------------------------------------------

cal_grid05a <- expand.grid(
  n = n_grid05a, x_type = x_grid05a, error_type = error_grid05a,
  null_class = null_class_grid05a,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

cal_grid05a$cell_id <- seq_len(nrow(cal_grid05a))

cal_grid05a$env_id <- as.integer(factor(paste(
  cal_grid05a$n, cal_grid05a$x_type, cal_grid05a$error_type, sep = "|")))


# ------------------------------------------------------------------------------
# Evaluation: 96 environments x 8 scenarios x 10 severities.
#
# seed_group_id EXCLUDES severity, so every severity level in an environment
# receives the same seed_base and CRN holds across the severity axis.
# ------------------------------------------------------------------------------

eval_rows <- list()

for (sc in scenario_grid05a) {
  sv <- severity_grids05a[[severity_family05a(sc)]]
  eval_rows[[length(eval_rows) + 1L]] <- expand.grid(
    n = n_grid05a, x_type = x_grid05a, error_type = error_grid05a,
    scenario = sc, severity_index = seq_along(sv),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
}

eval_grid05a <- do.call(rbind, eval_rows)

eval_grid05a$severity_target <- mapply(
  function(sc, ix) severity_grids05a[[severity_family05a(sc)]][ix],
  eval_grid05a$scenario, eval_grid05a$severity_index)

eval_grid05a$null_class <- vapply(eval_grid05a$scenario, null_class05a,
                                  character(1))

eval_grid05a$run_correct <-
  eval_grid05a$scenario %in% has_correct05a &
  eval_grid05a$severity_index %in% correct_state_idx05a

eval_grid05a$seed_group_id <- as.integer(factor(paste(
  eval_grid05a$n, eval_grid05a$x_type, eval_grid05a$error_type,
  eval_grid05a$scenario, sep = "|")))

eval_grid05a$env_id <- as.integer(factor(paste(
  eval_grid05a$n, eval_grid05a$x_type, eval_grid05a$error_type, sep = "|")))

eval_grid05a <- eval_grid05a[order(eval_grid05a$scenario,
                                   eval_grid05a$n,
                                   eval_grid05a$x_type,
                                   eval_grid05a$error_type,
                                   eval_grid05a$severity_index), ]

eval_grid05a$cell_id <- seq_len(nrow(eval_grid05a))

row.names(eval_grid05a) <- NULL


message("Calibration cells: ", nrow(cal_grid05a),
        "  (", sim_params05a$n_iters_cal, " draws each)")
message("Evaluation cells:  ", nrow(eval_grid05a),
        "  (", sim_params05a$n_iters_eval, " draws each)")
message("Environments:      ", length(unique(eval_grid05a$env_id)))
message("CRN seed groups:   ", length(unique(eval_grid05a$seed_group_id)))


# ==============================================================================
# 7. Directories
# ==============================================================================

output_dir <- "../output"

checkpoint_dir <- file.path(output_dir, "temp", "05a_detector")

for (dd in c(output_dir, file.path(output_dir, "temp"), checkpoint_dir)) {
  if (!dir.exists(dd)) dir.create(dd, recursive = TRUE, showWarnings = FALSE)
}

cutoff_path <- file.path(output_dir, "05a_null_cutoffs.rds")
grid_path   <- file.path(output_dir, "05a_design_grids.rds")

safe_save_rds(
  list(cal_grid = cal_grid05a, eval_grid = eval_grid05a,
       severity_grids = severity_grids05a, k_fractions = k_fraction_grid05a,
       params = sim_params05a, c0_pop = c0_pop05a,
       created = Sys.time()),
  grid_path)


# ==============================================================================
# 8. Workers
#
# Same architecture as 05: an explicit PSOCK cluster, C++ compiled once inside
# every worker, then future::plan(future::cluster).
# ==============================================================================

n_workers <- sim_params05a$workers

message("Starting ", n_workers, " parallel workers.")

cl <- parallel::makePSOCKcluster(n_workers)


# ------------------------------------------------------------------------------
# Frozen modules inside each worker
# ------------------------------------------------------------------------------

r_dir <- normalizePath("../R", mustWork = TRUE)

parallel::clusterExport(cl, "r_dir", envir = environment())

invisible(parallel::clusterCall(cl, function() {
  source(file.path(r_dir, "helpers_local.R"))
  source(file.path(r_dir, "diagnostics_classical.R"))
  source(file.path(r_dir, "dinkelbach_topk.R"))
  source(file.path(r_dir, "dgp_misspec_factory.R"))
  source(file.path(r_dir, "mis_sensitivity.R"))
  suppressPackageStartupMessages({
    requireNamespace("lmtest", quietly = TRUE)
    requireNamespace("sandwich", quietly = TRUE)
  })
  TRUE
}))


# ------------------------------------------------------------------------------
# 05a-local functions inside each worker
# ------------------------------------------------------------------------------

parallel::clusterExport(
  cl,
  c("calibrate_sd05a", "generate_misspec_data05a", "spec_tests05a",
    "fit_stats05a", "run_cal_iter05a", "run_eval_iter05a"),
  envir = environment())


# ------------------------------------------------------------------------------
# C++ MIS kernel
# ------------------------------------------------------------------------------

cpp_ok <- FALSE

if (isTRUE(sim_params05a$use_cpp) &&
    requireNamespace("Rcpp", quietly = TRUE) &&
    file.exists("../src/dinkelbach_topk_cpp.cpp")) {
  
  cpp_path <- normalizePath("../src/dinkelbach_topk_cpp.cpp", mustWork = TRUE)
  
  parallel::clusterExport(cl, "cpp_path", envir = environment())
  
  cpp_status <- parallel::clusterCall(cl, function() {
    if (!requireNamespace("Rcpp", quietly = TRUE)) return(FALSE)
    tryCatch({
      Rcpp::sourceCpp(cpp_path, rebuild = FALSE, showOutput = FALSE,
                      verbose = FALSE)
      exists("dinkelbach_topk_cpp", mode = "function")
    }, error = function(e) FALSE)
  })
  
  cpp_ok <- all(unlist(cpp_status))
}

message("C++ MIS kernel: ",
        if (cpp_ok) "enabled on all workers" else
          "disabled; using pure-R Dinkelbach")

future::plan(future::cluster, workers = cl)

message("Active workers: ", future::nbrOfWorkers())


# ==============================================================================
# 9. PHASE 1 -- CALIBRATION
#
# Disjoint seed stream. These draws never enter the evaluation sample, and the
# evaluation sample never contributes to a cutoff.
# ==============================================================================

message("\n", strrep("=", 78))
message("PHASE 1: CALIBRATION")
message(strrep("=", 78))

tryCatch({
  
  for (i in seq_len(nrow(cal_grid05a))) {
    
    cell <- cal_grid05a[i, , drop = FALSE]
    
    checkpoint_file <- file.path(
      checkpoint_dir, sprintf("cal_cell_%06d.rds", cell$cell_id))
    
    if (is_computed(checkpoint_file)) {
      message(sprintf("[cal %d/%d] cell %d done -- skipping.",
                      i, nrow(cal_grid05a), cell$cell_id))
      next
    }
    
    message(sprintf(
      "[cal %d/%d] N=%d | X=%s | error=%s | null=%s",
      i, nrow(cal_grid05a), cell$n, cell$x_type, cell$error_type,
      cell$null_class))
    
    seed_base <- as.integer(sim_params05a$seed_cal + 100000L * cell$cell_id)
    
    cell_results <- furrr::future_map_dfr(
      seq_len(sim_params05a$n_iters_cal),
      function(it) {
        run_cal_iter05a(
          iter = it, n = cell$n, x_type = cell$x_type,
          error_type = cell$error_type, null_class = cell$null_class,
          k_fractions = k_fraction_grid05a, seed_base = seed_base,
          use_cpp = cpp_ok, params = sim_params05a,
          c0 = c0_pop05a[[cell$x_type]])
      },
      .options = furrr::furrr_options(seed = TRUE, scheduling = 2))
    
    if (nrow(cell_results) == 0L) {
      warning("Calibration cell ", cell$cell_id, " produced no rows.")
      next
    }
    
    cell_results$cell_id    <- cell$cell_id
    cell_results$env_id     <- cell$env_id
    cell_results$n          <- cell$n
    cell_results$x_type     <- cell$x_type
    cell_results$error_type <- cell$error_type
    cell_results$null_class <- cell$null_class
    cell_results$phase      <- "calibration"
    
    safe_save_rds(cell_results, checkpoint_file)
    
    message(sprintf("  cal cell %d saved (%d rows).",
                    cell$cell_id, nrow(cell_results)))
  }
  
}, error = function(e) {
  message("CALIBRATION ERROR: ", conditionMessage(e))
  try(parallel::stopCluster(cl), silent = TRUE)
  stop(e)
})


# ==============================================================================
# 10. NULL CUTOFFS
#
# Influence statistics: 95th percentile.
# Specification tests: 5th percentile of the null p-value distribution, which
#   is the size correction. Pilots showed RESET running at nominal size
#   0.085-0.110, so comparing it against MIS at nominal 0.05 would be unfair.
# ==============================================================================

message("\n", strrep("=", 78))
message("COMPUTING NULL CUTOFFS")
message(strrep("=", 78))

cut_rows <- list()

for (i in seq_len(nrow(cal_grid05a))) {
  
  f <- file.path(checkpoint_dir,
                 sprintf("cal_cell_%06d.rds", cal_grid05a$cell_id[i]))
  
  if (!is_computed(f)) {
    warning("Missing calibration checkpoint: ", f)
    next
  }
  
  dd <- readRDS(f)
  
  for (kf in k_fraction_grid05a) {
    
    s <- dd[dd$k_fraction == kf, , drop = FALSE]
    if (nrow(s) == 0L) next
    
    cut_rows[[length(cut_rows) + 1L]] <- data.frame(
      n = s$n[1], x_type = s$x_type[1], error_type = s$error_type[1],
      null_class = s$null_class[1], k_fraction = kf,
      n_cal = nrow(s),
      cut_MIS      = stats::quantile(s$MIS,      0.95, na.rm = TRUE, names = FALSE),
      cut_Cook     = stats::quantile(s$Cook,     0.95, na.rm = TRUE, names = FALSE),
      cut_DFBETAS  = stats::quantile(s$DFBETAS,  0.95, na.rm = TRUE, names = FALSE),
      cut_Leverage = stats::quantile(s$Leverage, 0.95, na.rm = TRUE, names = FALSE),
      cut_reset_p  = stats::quantile(s$reset_p,  0.05, na.rm = TRUE, names = FALSE),
      cut_bp_p     = stats::quantile(s$bp_p,     0.05, na.rm = TRUE, names = FALSE),
      null_mean_MIS = mean(s$MIS, na.rm = TRUE),
      finite_MIS_rate = mean(is.finite(s$MIS)),
      stringsAsFactors = FALSE)
  }
}

null_cutoffs05a <- do.call(rbind, cut_rows)

safe_save_rds(null_cutoffs05a, cutoff_path)

message("Cutoffs saved: ", nrow(null_cutoffs05a), " rows -> ", cutoff_path)

# Quick sanity read on the sqrt(n) property established in pilot 3.
if (!is.null(null_cutoffs05a)) {
  chk <- null_cutoffs05a[null_cutoffs05a$k_fraction == 0.05 &
                           null_cutoffs05a$null_class == "x" &
                           null_cutoffs05a$x_type == "normal" &
                           null_cutoffs05a$error_type == "normal", ]
  if (nrow(chk) > 1L) {
    sl <- stats::coef(stats::lm(log(null_mean_MIS) ~ log(n), data = chk))[2]
    message("Null MIS log-log slope on n (expect ~0.5): ", round(sl, 4))
  }
}


# ==============================================================================
# 11. PHASE 2 -- EVALUATION
#
# CRN across severity is retained here: seed_group_id excludes severity, so
# every severity level in an environment reproduces the same X, eps0, z and g.
# ==============================================================================

message("\n", strrep("=", 78))
message("PHASE 2: EVALUATION")
message(strrep("=", 78))

tryCatch({
  
  for (i in seq_len(nrow(eval_grid05a))) {
    
    cell <- eval_grid05a[i, , drop = FALSE]
    
    checkpoint_file <- file.path(
      checkpoint_dir, sprintf("eval_cell_%06d.rds", cell$cell_id))
    
    if (is_computed(checkpoint_file)) {
      if (i %% 200L == 0L)
        message(sprintf("[eval %d/%d] cell %d done -- skipping.",
                        i, nrow(eval_grid05a), cell$cell_id))
      next
    }
    
    message(sprintf(
      paste0("[eval %d/%d] N=%d | X=%s | error=%s | scenario=%s | ",
             "sev=%d (%.3f) | seed-group=%d%s"),
      i, nrow(eval_grid05a), cell$n, cell$x_type, cell$error_type,
      cell$scenario, cell$severity_index, cell$severity_target,
      cell$seed_group_id,
      if (isTRUE(cell$run_correct)) " | +correct" else ""))
    
    seed_base <- as.integer(sim_params05a$seed_eval +
                              100000L * cell$seed_group_id)
    
    cell_results <- furrr::future_map_dfr(
      seq_len(sim_params05a$n_iters_eval),
      function(it) {
        run_eval_iter05a(
          iter = it, n = cell$n, x_type = cell$x_type,
          error_type = cell$error_type, scenario = cell$scenario,
          severity_index = cell$severity_index,
          severity_target = cell$severity_target,
          k_fractions = k_fraction_grid05a, seed_base = seed_base,
          use_cpp = cpp_ok, params = sim_params05a,
          c0 = c0_pop05a[[cell$x_type]],
          run_correct = cell$run_correct)
      },
      .options = furrr::furrr_options(seed = TRUE, scheduling = 2))
    
    if (nrow(cell_results) == 0L) {
      warning("Evaluation cell ", cell$cell_id, " produced no rows.")
      next
    }
    
    cell_results$cell_id         <- cell$cell_id
    cell_results$seed_group_id   <- cell$seed_group_id
    cell_results$env_id          <- cell$env_id
    cell_results$n               <- cell$n
    cell_results$x_type          <- cell$x_type
    cell_results$error_type      <- cell$error_type
    cell_results$scenario        <- cell$scenario
    cell_results$severity_index  <- cell$severity_index
    cell_results$severity_target <- cell$severity_target
    cell_results$null_class      <- cell$null_class
    cell_results$phase           <- "evaluation"
    
    safe_save_rds(cell_results, checkpoint_file)
    
    if (i %% 25L == 0L || i == nrow(eval_grid05a))
      message(sprintf("  eval cell %d saved (%d rows). [%.1f%% complete]",
                      cell$cell_id, nrow(cell_results),
                      100 * i / nrow(eval_grid05a)))
  }
  
}, error = function(e) {
  message("EVALUATION ERROR: ", conditionMessage(e))
  try(parallel::stopCluster(cl), silent = TRUE)
  stop(e)
})


# ==============================================================================
# 12. Shutdown
# ==============================================================================

future::plan(future::sequential)

try(parallel::stopCluster(cl), silent = TRUE)

n_cal_done  <- length(list.files(checkpoint_dir, pattern = "^cal_cell_"))
n_eval_done <- length(list.files(checkpoint_dir, pattern = "^eval_cell_"))

message("\n", strrep("=", 78))
message("05a COMPLETE")
message(strrep("=", 78))
message("Calibration checkpoints: ", n_cal_done, " / ", nrow(cal_grid05a))
message("Evaluation checkpoints:  ", n_eval_done, " / ", nrow(eval_grid05a))
message("Cutoffs:                 ", cutoff_path)
message("Design grids:            ", grid_path)
message("Checkpoint directory:    ", checkpoint_dir)
message("")
message("Frozen 05 output was not read, written, or modified.")
message("Next: 05b post-processing (size, power, isotonic 80% boundary).")
message(strrep("=", 78))