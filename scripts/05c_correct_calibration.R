# ==============================================================================
# File: /scripts/05c_correct_calibration.R
#
# CORRECT-MODEL NULL CALIBRATION
#
# WHY THIS EXISTS
#
# Every cutoff in 05a_null_cutoffs.rds was estimated under a WRONG-model null
# class: y ~ x, or y ~ x + z for missing_interaction. The correct-state rows in
# 05a were then scored against those same cutoffs. But the correct models carry
# extra regressors, and the null level of standardized delta rises with the
# number of regressors -- so the reported false-positive rate after
# respecification was inflated wherever the correct model is much larger than
# the wrong one:
#
#   scenario              power_correct (05b)   correct formula
#   nonlinear                     0.0387        y ~ x + I(x^2)
#   missing_interaction           0.0643        y ~ x + z + x:z
#   ovb                           0.0188        y ~ x + z
#   heterogeneous                 0.3979        y ~ x + g + x:g     <- inflated
#   threshold                     0.3716        y ~ x + hinge       <- inflated
#
# This script estimates the null distribution UNDER EACH CORRECT MODEL and
# recomputes the specificity table against the right cutoffs.
#
# mean_reduction is cutoff-free and will not move. The threshold anomaly
# (-7.97 under contaminated X) is therefore expected to survive; what changes
# is whether the calibrated TEST still holds its size there, which is a much
# weaker and more reportable failure than an inflated rejection rate.
#
# ------------------------------------------------------------------------------
# DESIGN
#
# At severity 0 the DGP is identical for every scenario, so ONE null draw per
# (n, x_type, error_type) serves all five correct models. The draw is generated
# through scenario = "threshold" at sev = 0 -- not "correct" -- because that is
# the only path that populates the hinge column, and z and g are generated
# unconditionally on every path. RNG consumption therefore matches the
# evaluation draws exactly.
#
#   96 environments x 2500 draws, five formulas fitted per draw.
#
# Classical diagnostics and specification tests are off by default: the
# specificity exhibit is about MIS, and skipping them cuts roughly 60% of the
# cost. Set do_classical05c <- TRUE if you want the full set.
#
# COST: ~2.5 h on 12 workers with classical off; ~6 h with it on.
#
# OUTPUT
#   ../output/temp/05c_correct_cal/cal_cell_XXXXXX.rds   (96 files)
#   ../output/05a_correct_cutoffs.rds
#   ../output/05a_specificity_corrected.rds
#
# Frozen 05 and the 05a checkpoints are read-only here.
# ==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(future)
  library(furrr)
  library(parallel)
})

source("../R/helpers_local.R")
source("../R/utils_checkpoint.R")
source("../R/diagnostics_classical.R")
source("../R/dinkelbach_topk.R")
source("../R/dgp_misspec_factory.R")
source("../R/mis_sensitivity.R")

output_dir     <- "../output"
cal_dir        <- file.path(output_dir, "temp", "05c_correct_cal")
eval_ckpt_dir  <- file.path(output_dir, "temp", "05a_detector")

if (!dir.exists(cal_dir)) dir.create(cal_dir, recursive = TRUE,
                                     showWarnings = FALSE)

grids <- readRDS(file.path(output_dir, "05a_design_grids.rds"))

k_grid    <- grids$k_fractions
params    <- grids$params
c0_pop    <- grids$c0_pop
eval_grid <- grids$eval_grid

n_iters05c    <- 2500L
seed_cal05c   <- 88260828L      # disjoint from 05a's 20260828 / 77260828
do_classical05c <- FALSE
PRIMARY_K     <- 0.05


# ==============================================================================
# 1. The five correct-model classes
# ==============================================================================

correct_formulas05c <- list(
  nonlinear           = y ~ x + I(x^2),
  threshold           = y ~ x + hinge,
  ovb                 = y ~ x + z,
  heterogeneous       = y ~ x + g + x:g,
  missing_interaction = y ~ x + z + x:z
)


# ==============================================================================
# 2. DGP -- byte-identical to the 05a definition
#
# Reproduced here rather than sourced so that 05a's script is not made a
# dependency of this one. If 05a's generator is ever edited, edit both.
# ==============================================================================

calibrate_sd05a <- function(component, error, target_ratio) {
  if (!is.finite(target_ratio) || target_ratio <= 0)
    return(list(term = rep(0, length(component)), raw_coef = 0,
                realized_ratio = 0))
  sc_comp <- stats::sd(component, na.rm = TRUE)
  sc_err  <- stats::sd(error, na.rm = TRUE)
  if (!is.finite(sc_comp) || sc_comp <= .Machine$double.eps)
    return(list(term = rep(0, length(component)), raw_coef = 0,
                realized_ratio = 0))
  raw_coef <- target_ratio * sc_err / sc_comp
  term     <- raw_coef * component
  list(term = term, raw_coef = raw_coef,
       realized_ratio = stats::sd(term, na.rm = TRUE) / sc_err)
}

generate_null_draw05c <- function(n, x_type, error_type, c0_pop,
                                  mix_prop = 0.10, rho_xz = 0.50,
                                  hetero_group_prop = 0.25,
                                  gpd_shape = 0.25, pareto_shape = 3,
                                  skew_t_df = 5) {
  
  x <- generate_vector05(n, x_type, mix_prop = mix_prop,
                         skew_t_df = skew_t_df, gpd_shape = gpd_shape,
                         pareto_shape = pareto_shape, center = FALSE)
  
  eps0 <- generate_vector05(n, error_type, mix_prop = mix_prop,
                            skew_t_df = skew_t_df, gpd_shape = gpd_shape,
                            pareto_shape = pareto_shape, center = TRUE)
  
  # Same order and count of RNG calls as the evaluation generator.
  z <- generate_correlated_z05(x, rho = rho_xz)
  g <- as.integer(stats::runif(n) < hetero_group_prop)
  
  hinge <- pmax(x - c0_pop, 0)
  
  data.frame(y = as.numeric(x + eps0), x = as.numeric(x),
             z = as.numeric(z), g = g, hinge = as.numeric(hinge))
}


# ==============================================================================
# 3. Statistics for one fitted correct model
# ==============================================================================

fit_stats05c <- function(formula, dat, k_fractions, use_cpp, do_classical) {
  
  n <- nrow(dat)
  
  mod <- tryCatch(stats::lm(formula, data = dat), error = function(e) NULL)
  if (is.null(mod)) return(NULL)
  
  sm <- tryCatch(summary(mod)$coefficients, error = function(e) NULL)
  if (is.null(sm) || !("x" %in% rownames(sm))) return(NULL)
  
  b <- unname(sm["x", "Estimate"]); s <- unname(sm["x", "Std. Error"])
  if (!is.finite(b) || !is.finite(s) || s <= 0) return(NULL)
  
  target_pos <- match("x", rownames(sm))
  
  rk <- NULL
  if (isTRUE(do_classical)) {
    cl <- tryCatch(get_all_classical(mod, "x"), error = function(e) NULL)
    if (!is.null(cl)) rk <- list(
      Cook     = cl$id[order(cl$cooks_d, decreasing = TRUE)],
      DFBETAS  = cl$id[order(abs(cl$dfbetas_target), decreasing = TRUE)],
      Leverage = cl$id[order(cl$leverage, decreasing = TRUE)])
  }
  
  rows <- vector("list", length(k_fractions))
  
  for (j in seq_along(k_fractions)) {
    
    kf <- k_fractions[j]
    k  <- max(1L, as.integer(round(kf * n)))
    
    m <- tryCatch(
      run_mis_sensitivity(mod_full = mod, formula = formula, data = dat,
                          k = k, target_var = "x", target_pos = target_pos,
                          use_cpp = use_cpp),
      error = function(e) NULL)
    
    cls <- c(Cook = NA_real_, DFBETAS = NA_real_, Leverage = NA_real_)
    if (!is.null(rk)) {
      for (nm in names(rk)) {
        f <- tryCatch(fit_target_ols05(formula, dat, "x", rk[[nm]][seq_len(k)]),
                      error = function(e) NULL)
        if (!is.null(f)) cls[nm] <- unname(abs(f["coef"] - b) / s)
      }
    }
    
    rows[[j]] <- data.frame(
      k_fraction = kf, k = k, full_coef = b, full_se = s,
      MIS = if (is.null(m)) NA_real_ else m$standardized_delta,
      Cook = unname(cls["Cook"]), DFBETAS = unname(cls["DFBETAS"]),
      Leverage = unname(cls["Leverage"]),
      stringsAsFactors = FALSE)
  }
  
  do.call(rbind, rows)
}


# ==============================================================================
# 4. One calibration draw -- all five correct models on the same null data
# ==============================================================================

run_cal_iter05c <- function(iter, n, x_type, error_type, k_fractions,
                            seed_base, use_cpp, params, c0, formulas,
                            do_classical) {
  
  set.seed(seed_base + iter)
  
  dat <- tryCatch(
    generate_null_draw05c(n = n, x_type = x_type, error_type = error_type,
                          c0_pop = c0, mix_prop = params$mix_prop,
                          rho_xz = params$rho_xz,
                          hetero_group_prop = params$hetero_group_prop,
                          gpd_shape = params$gpd_shape,
                          pareto_shape = params$pareto_shape,
                          skew_t_df = params$skew_t_df),
    error = function(e) NULL)
  
  if (is.null(dat)) return(NULL)
  
  out <- list()
  
  for (nm in names(formulas)) {
    st <- fit_stats05c(formulas[[nm]], dat, k_fractions, use_cpp, do_classical)
    if (!is.null(st)) {
      st$correct_class <- nm
      out[[length(out) + 1L]] <- st
    }
  }
  
  if (length(out) == 0L) return(NULL)
  
  res <- do.call(rbind, out)
  res$iter <- iter
  res
}


# ==============================================================================
# 5. Calibration grid -- 96 environments, no scenario axis
# ==============================================================================

cal_grid05c <- unique(eval_grid[, c("n", "x_type", "error_type", "env_id")])
cal_grid05c <- cal_grid05c[order(cal_grid05c$env_id), ]
cal_grid05c$cell_id <- seq_len(nrow(cal_grid05c))
row.names(cal_grid05c) <- NULL

message("Correct-model calibration cells: ", nrow(cal_grid05c),
        " x ", n_iters05c, " draws x ", length(correct_formulas05c),
        " formulas")


# ==============================================================================
# 6. Workers -- same architecture as 05a
# ==============================================================================

n_workers <- params$workers
message("Starting ", n_workers, " workers.")
cl <- parallel::makePSOCKcluster(n_workers)

r_dir <- normalizePath("../R", mustWork = TRUE)
parallel::clusterExport(cl, "r_dir", envir = environment())

invisible(parallel::clusterCall(cl, function() {
  source(file.path(r_dir, "helpers_local.R"))
  source(file.path(r_dir, "diagnostics_classical.R"))
  source(file.path(r_dir, "dinkelbach_topk.R"))
  source(file.path(r_dir, "dgp_misspec_factory.R"))
  source(file.path(r_dir, "mis_sensitivity.R"))
  TRUE
}))

parallel::clusterExport(
  cl, c("calibrate_sd05a", "generate_null_draw05c", "fit_stats05c",
        "run_cal_iter05c"), envir = environment())

cpp_ok <- FALSE
if (isTRUE(params$use_cpp) && requireNamespace("Rcpp", quietly = TRUE) &&
    file.exists("../src/dinkelbach_topk_cpp.cpp")) {
  cpp_path <- normalizePath("../src/dinkelbach_topk_cpp.cpp", mustWork = TRUE)
  parallel::clusterExport(cl, "cpp_path", envir = environment())
  cpp_ok <- all(unlist(parallel::clusterCall(cl, function() {
    if (!requireNamespace("Rcpp", quietly = TRUE)) return(FALSE)
    tryCatch({
      Rcpp::sourceCpp(cpp_path, rebuild = FALSE, showOutput = FALSE,
                      verbose = FALSE)
      exists("dinkelbach_topk_cpp", mode = "function")
    }, error = function(e) FALSE)
  })))
}

message("C++ MIS kernel: ", if (cpp_ok) "enabled" else "disabled")

future::plan(future::cluster, workers = cl)


# ==============================================================================
# 7. Run
# ==============================================================================

tryCatch({
  
  for (i in seq_len(nrow(cal_grid05c))) {
    
    cell <- cal_grid05c[i, , drop = FALSE]
    f <- file.path(cal_dir, sprintf("cal_cell_%06d.rds", cell$cell_id))
    
    if (is_computed(f)) {
      message(sprintf("[%d/%d] cell %d done -- skipping.",
                      i, nrow(cal_grid05c), cell$cell_id))
      next
    }
    
    message(sprintf("[%d/%d] N=%d | X=%s | error=%s",
                    i, nrow(cal_grid05c), cell$n, cell$x_type,
                    cell$error_type))
    
    seed_base <- as.integer(seed_cal05c + 100000L * cell$cell_id)
    
    res <- furrr::future_map_dfr(
      seq_len(n_iters05c),
      function(it) run_cal_iter05c(
        iter = it, n = cell$n, x_type = cell$x_type,
        error_type = cell$error_type, k_fractions = k_grid,
        seed_base = seed_base, use_cpp = cpp_ok, params = params,
        c0 = c0_pop[[cell$x_type]], formulas = correct_formulas05c,
        do_classical = do_classical05c),
      .options = furrr::furrr_options(seed = TRUE, scheduling = 2))
    
    if (nrow(res) == 0L) { warning("cell ", cell$cell_id, " empty."); next }
    
    res$cell_id    <- cell$cell_id
    res$env_id     <- cell$env_id
    res$n          <- cell$n
    res$x_type     <- cell$x_type
    res$error_type <- cell$error_type
    
    safe_save_rds(res, f)
    message(sprintf("  saved (%d rows).", nrow(res)))
  }
  
}, error = function(e) {
  try(parallel::stopCluster(cl), silent = TRUE)
  stop(e)
})

future::plan(future::sequential)
try(parallel::stopCluster(cl), silent = TRUE)


# ==============================================================================
# 8. Correct-model cutoffs
# ==============================================================================

message("\nComputing correct-model cutoffs ...")

cut_rows <- list()

for (i in seq_len(nrow(cal_grid05c))) {
  
  f <- file.path(cal_dir, sprintf("cal_cell_%06d.rds", cal_grid05c$cell_id[i]))
  if (!is_computed(f)) { warning("missing ", f); next }
  
  dd <- readRDS(f)
  
  for (cc in unique(dd$correct_class)) {
    for (kf in k_grid) {
      s <- dd[dd$correct_class == cc & dd$k_fraction == kf, , drop = FALSE]
      if (nrow(s) == 0L) next
      cut_rows[[length(cut_rows) + 1L]] <- data.frame(
        n = s$n[1], x_type = s$x_type[1], error_type = s$error_type[1],
        correct_class = cc, k_fraction = kf, n_cal = nrow(s),
        cut_MIS = stats::quantile(s$MIS, 0.95, na.rm = TRUE, names = FALSE),
        null_mean_MIS = mean(s$MIS, na.rm = TRUE),
        finite_MIS_rate = mean(is.finite(s$MIS)),
        stringsAsFactors = FALSE)
    }
  }
}

correct_cutoffs <- do.call(rbind, cut_rows)
safe_save_rds(correct_cutoffs,
              file.path(output_dir, "05a_correct_cutoffs.rds"))

cat("\n", strrep("=", 90), "\n", sep = "")
cat("NULL LEVEL: WRONG MODEL VERSUS CORRECT MODEL (k = 5%, mean over ",
    "environments)\n", sep = "")
cat(strrep("=", 90), "\n", sep = "")

wrong_cut <- readRDS(file.path(output_dir, "05a_null_cutoffs.rds"))
wc <- wrong_cut[wrong_cut$k_fraction == PRIMARY_K &
                  wrong_cut$null_class == "x", ]
cc5 <- correct_cutoffs[correct_cutoffs$k_fraction == PRIMARY_K, ]

cmp <- do.call(rbind, lapply(split(cc5, cc5$correct_class), function(g)
  data.frame(correct_class = g$correct_class[1],
             correct_null_mean = mean(g$null_mean_MIS, na.rm = TRUE),
             correct_cutoff = mean(g$cut_MIS, na.rm = TRUE),
             wrong_null_mean = mean(wc$null_mean_MIS, na.rm = TRUE),
             wrong_cutoff = mean(wc$cut_MIS, na.rm = TRUE),
             stringsAsFactors = FALSE)))
cmp$cutoff_ratio <- cmp$correct_cutoff / cmp$wrong_cutoff
print(cmp, digits = 5, row.names = FALSE)

cat("\ncutoff_ratio > 1 quantifies exactly how much the 05b specificity table\n")
cat("was over-rejecting by scoring correct-model draws against y ~ x cutoffs.\n")


# ==============================================================================
# 9. Corrected specificity table
#
# Re-streams only the correct-state rows from the 05a evaluation checkpoints
# (5 scenarios x 4 severity indices x 96 environments = 1920 cells).
# ==============================================================================

message("\nRecomputing specificity against correct-model cutoffs ...")

ck <- paste(correct_cutoffs$n, correct_cutoffs$x_type,
            correct_cutoffs$error_type, correct_cutoffs$correct_class,
            correct_cutoffs$k_fraction, sep = "|")

target_cells <- eval_grid[eval_grid$run_correct, ]

spec_rows <- list()

for (i in seq_len(nrow(target_cells))) {
  
  cell <- target_cells[i, , drop = FALSE]
  f <- file.path(eval_ckpt_dir, sprintf("eval_cell_%06d.rds", cell$cell_id))
  if (!is_computed(f)) next
  
  dd <- readRDS(f)
  cr <- dd[dd$model_state == "correct", , drop = FALSE]
  wr <- dd[dd$model_state == "wrong", , drop = FALSE]
  if (nrow(cr) == 0L) next
  
  for (kf in k_grid) {
    
    sc <- cr[cr$k_fraction == kf, , drop = FALSE]
    sw <- wr[wr$k_fraction == kf, , drop = FALSE]
    if (nrow(sc) == 0L) next
    
    j <- match(paste(cell$n, cell$x_type, cell$error_type,
                     cell$scenario, kf, sep = "|"), ck)
    if (is.na(j)) next
    
    spec_rows[[length(spec_rows) + 1L]] <- data.frame(
      cell_id = cell$cell_id, env_id = cell$env_id, n = cell$n,
      x_type = cell$x_type, error_type = cell$error_type,
      scenario = cell$scenario, severity_index = cell$severity_index,
      severity_target = cell$severity_target, k_fraction = kf,
      mean_MIS_correct = mean(sc$MIS, na.rm = TRUE),
      mean_MIS_wrong = if (nrow(sw)) mean(sw$MIS, na.rm = TRUE) else NA_real_,
      # scored against the CORRECT-model cutoff
      power_correct_calibrated =
        mean(sc$MIS > correct_cutoffs$cut_MIS[j], na.rm = TRUE),
      correct_cutoff = correct_cutoffs$cut_MIS[j],
      stringsAsFactors = FALSE)
  }
  
  if (i %% 250L == 0L)
    message(sprintf("  %d / %d", i, nrow(target_cells)))
}

spec <- do.call(rbind, spec_rows)
spec$reduction <- spec$mean_MIS_wrong - spec$mean_MIS_correct

safe_save_rds(spec, file.path(output_dir, "05a_specificity_corrected.rds"))

p5 <- spec[spec$k_fraction == PRIMARY_K, ]

cat("\n", strrep("=", 90), "\n", sep = "")
cat("CORRECTED SPECIFICITY -- k = 5%\n")
cat(strrep("=", 90), "\n", sep = "")

cat("\nSize under correct specification, severity 0 (should be ~0.05):\n\n")
s0 <- p5[p5$severity_index == 1L, ]
print(aggregate(list(size_correct = s0$power_correct_calibrated),
                by = list(scenario = s0$scenario), FUN = mean),
      digits = 5, row.names = FALSE)

cat("\nFalse-positive rate under correct specification, positive severities:\n\n")
pp <- p5[p5$severity_index > 1L, ]
print(aggregate(
  list(power_correct_calibrated = pp$power_correct_calibrated,
       mean_reduction = pp$reduction,
       prop_correct_lower = as.numeric(pp$reduction > 0)),
  by = list(scenario = pp$scenario), FUN = function(v) mean(v, na.rm = TRUE)),
  digits = 5, row.names = FALSE)

cat("\nCompare power_correct_calibrated against 05b's uncalibrated values:\n")
cat("  heterogeneous 0.3979 | threshold 0.3716 | nonlinear 0.0387\n")
cat("  missing_interaction 0.0643 | ovb 0.0188\n")

cat("\nThreshold by X distribution -- the anomaly localises here:\n\n")
th <- pp[pp$scenario == "threshold", ]
print(aggregate(
  list(power_correct = th$power_correct_calibrated,
       mean_reduction = th$reduction),
  by = list(x_type = th$x_type), FUN = mean), digits = 5, row.names = FALSE)

cat("\nIf power_correct is near 0.05 for threshold even where mean_reduction\n")
cat("is negative, then the statistic inflates under over-specification with\n")
cat("heavy-tailed X while the calibrated TEST still holds its size. That is\n")
cat("a far weaker caveat than an inflated rejection rate.\n")

cat("\nSaved: 05a_correct_cutoffs.rds, 05a_specificity_corrected.rds\n")
cat(strrep("=", 90), "\n", sep = "")