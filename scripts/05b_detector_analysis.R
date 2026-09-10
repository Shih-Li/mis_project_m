# ==============================================================================
# File: /scripts/05b_detector_analysis.R
#
# POST-PROCESSING FOR THE 05a DETECTOR EXPERIMENT
#
# Two stages:
#
#   STAGE 1  Stream the 7680 evaluation checkpoints one file at a time, join
#            each to its independently-estimated null cutoff, and reduce to a
#            compact per-(cell, k, model_state) summary. Saved to disk, so
#            stage 2 can be re-run instantly while drafting exhibits.
#
#   STAGE 2  Build the report: size, power, MIS vs classical influence
#            diagnostics, MIS vs size-corrected specification tests, isotonic
#            80% boundary, effect of k, respecification, localisation, the
#            endogeneity invariance result, and environment heterogeneity.
#
# MEMORY: the raw checkpoints hold roughly 20 million rows. Nothing here ever
# binds more than one file at a time.
#
# Reads only from ../output/temp/05a_detector and ../output/05a_*.rds.
# Frozen 05 output is not touched.
# ==============================================================================

rm(list = ls())

source("../R/utils_checkpoint.R")

output_dir     <- "../output"
checkpoint_dir <- file.path(output_dir, "temp", "05a_detector")

summary_path <- file.path(output_dir, "05a_detector_summary.rds")
bound_path   <- file.path(output_dir, "05a_detection_boundary.rds")

grids   <- readRDS(file.path(output_dir, "05a_design_grids.rds"))
cutoffs <- readRDS(file.path(output_dir, "05a_null_cutoffs.rds"))

eval_grid <- grids$eval_grid
k_grid    <- grids$k_fractions
sev_grids <- grids$severity_grids
alpha     <- grids$params$alpha

PRIMARY_K <- 0.05

line   <- strrep("=", 100)
subline <- strrep("-", 100)

pr <- function(x, digits = 6) print(as.data.frame(x), digits = digits,
                                    row.names = FALSE)


# ==============================================================================
# STAGE 1 -- reduce checkpoints to a per-cell summary
# ==============================================================================

if (!file.exists(summary_path)) {
  
  message("STAGE 1: reducing ", nrow(eval_grid), " checkpoints ...")
  
  cut_key <- paste(cutoffs$n, cutoffs$x_type, cutoffs$error_type,
                   cutoffs$null_class, cutoffs$k_fraction, sep = "|")
  
  rows <- vector("list", nrow(eval_grid) * length(k_grid) * 2L)
  ptr  <- 0L
  missing_files <- integer(0)
  
  for (i in seq_len(nrow(eval_grid))) {
    
    cell <- eval_grid[i, , drop = FALSE]
    f <- file.path(checkpoint_dir, sprintf("eval_cell_%06d.rds", cell$cell_id))
    
    if (!is_computed(f)) { missing_files <- c(missing_files, cell$cell_id); next }
    
    dd <- readRDS(f)
    
    for (st in unique(dd$model_state)) {
      for (kf in k_grid) {
        
        s <- dd[dd$model_state == st & dd$k_fraction == kf, , drop = FALSE]
        if (nrow(s) == 0L) next
        
        ci <- match(paste(cell$n, cell$x_type, cell$error_type,
                          cell$null_class, kf, sep = "|"), cut_key)
        if (is.na(ci)) next
        
        cMIS <- cutoffs$cut_MIS[ci]
        cCk  <- cutoffs$cut_Cook[ci]
        cDF  <- cutoffs$cut_DFBETAS[ci]
        cLv  <- cutoffs$cut_Leverage[ci]
        cRS  <- cutoffs$cut_reset_p[ci]
        cBP  <- cutoffs$cut_bp_p[ci]
        
        rMIS <- s$MIS      > cMIS
        rCk  <- s$Cook     > cCk
        rDF  <- s$DFBETAS  > cDF
        rLv  <- s$Leverage > cLv
        rRS  <- s$reset_p  < cRS
        rBP  <- s$bp_p     < cBP
        
        ptr <- ptr + 1L
        rows[[ptr]] <- data.frame(
          cell_id = cell$cell_id, env_id = cell$env_id,
          seed_group_id = cell$seed_group_id,
          n = cell$n, x_type = cell$x_type, error_type = cell$error_type,
          scenario = cell$scenario, severity_index = cell$severity_index,
          severity_target = cell$severity_target,
          model_state = st, k_fraction = kf,
          draws = nrow(s),
          
          finite_MIS = mean(is.finite(s$MIS)),
          
          mean_MIS = mean(s$MIS, na.rm = TRUE),
          median_MIS = stats::median(s$MIS, na.rm = TRUE),
          mean_Cook = mean(s$Cook, na.rm = TRUE),
          mean_DFBETAS = mean(s$DFBETAS, na.rm = TRUE),
          mean_Leverage = mean(s$Leverage, na.rm = TRUE),
          
          rej_MIS = mean(rMIS, na.rm = TRUE),
          rej_Cook = mean(rCk, na.rm = TRUE),
          rej_DFBETAS = mean(rDF, na.rm = TRUE),
          rej_Leverage = mean(rLv, na.rm = TRUE),
          
          rej_RESET_cal = mean(s$reset_p < cRS, na.rm = TRUE),
          rej_BP_cal = mean(s$bp_p < cBP, na.rm = TRUE),
          rej_RESET_nom = mean(s$reset_p < alpha, na.rm = TRUE),
          rej_BP_nom = mean(s$bp_p < alpha, na.rm = TRUE),
          
          # McNemar-style net advantage: MIS alone minus rival alone
          net_vs_Cook = mean(rMIS & !rCk, na.rm = TRUE) -
            mean(!rMIS & rCk, na.rm = TRUE),
          net_vs_DFBETAS = mean(rMIS & !rDF, na.rm = TRUE) -
            mean(!rMIS & rDF, na.rm = TRUE),
          net_vs_Leverage = mean(rMIS & !rLv, na.rm = TRUE) -
            mean(!rMIS & rLv, na.rm = TRUE),
          net_vs_RESET = mean(rMIS & !(s$reset_p < cRS), na.rm = TRUE) -
            mean(!rMIS & (s$reset_p < cRS), na.rm = TRUE),
          net_vs_BP = mean(rMIS & !(s$bp_p < cBP), na.rm = TRUE) -
            mean(!rMIS & (s$bp_p < cBP), na.rm = TRUE),
          # Components of the net advantage. net_* is a difference and cannot
          # distinguish "MIS supersedes the rival" from "MIS complements it":
          # +0.30 is 0.35/0.05 or 0.60/0.30, which are different claims.
          only_MIS_vs_Cook     = mean(rMIS & !rCk, na.rm = TRUE),
          only_Cook_vs_MIS     = mean(!rMIS & rCk, na.rm = TRUE),
          only_MIS_vs_DFBETAS  = mean(rMIS & !rDF, na.rm = TRUE),
          only_DFBETAS_vs_MIS  = mean(!rMIS & rDF, na.rm = TRUE),
          only_MIS_vs_Leverage = mean(rMIS & !rLv, na.rm = TRUE),
          only_Leverage_vs_MIS = mean(!rMIS & rLv, na.rm = TRUE),
          only_MIS_vs_RESET    = mean(rMIS & !rRS, na.rm = TRUE),
          only_RESET_vs_MIS    = mean(!rMIS & rRS, na.rm = TRUE),
          only_MIS_vs_BP       = mean(rMIS & !rBP, na.rm = TRUE),
          only_BP_vs_MIS       = mean(!rMIS & rBP, na.rm = TRUE),
          mean_precision = mean(s$mis_precision, na.rm = TRUE),
          mean_recall = mean(s$mis_recall, na.rm = TRUE),
          mean_lift = mean(s$mis_lift, na.rm = TRUE),
          mean_affected = mean(s$affected_fraction, na.rm = TRUE),
          
          # CRN probe: under linear endogeneity the statistic is invariant to
          # severity, so its within-cell SD across draws is the only variation.
          sd_MIS = stats::sd(s$MIS, na.rm = TRUE),
          
          stringsAsFactors = FALSE)
      }
    }
    
    if (i %% 500L == 0L)
      message(sprintf("  %d / %d (%.1f%%)", i, nrow(eval_grid),
                      100 * i / nrow(eval_grid)))
  }
  
  det_summary <- do.call(rbind, rows[seq_len(ptr)])
  attr(det_summary, "missing_cells") <- missing_files
  
  safe_save_rds(det_summary, summary_path)
  message("Summary saved: ", nrow(det_summary), " rows -> ", summary_path)
  
} else {
  message("Loading existing summary ...")
  det_summary <- readRDS(summary_path)
}

S <- det_summary
W <- S[S$model_state == "wrong", ]
P <- W[W$k_fraction == PRIMARY_K, ]


# ==============================================================================
# Helper: isotonic 80% boundary
#
# Raw power curves are non-monotone at finite draws, so fit an isotonic curve
# before interpolating the first crossing. Never crossing is right-censored.
# ==============================================================================

boundary80 <- function(sev, pow, target = 0.80) {
  ok <- is.finite(pow) & is.finite(sev)
  if (sum(ok) < 2L) return(NA_real_)
  s <- sev[ok]; p <- pow[ok]; o <- order(s); s <- s[o]; p <- p[o]
  yf <- stats::isoreg(s, p)$yf
  if (all(yf < target)) return(NA_real_)
  if (yf[1] >= target) return(s[1])
  j <- which(yf >= target)[1]
  if (yf[j] == yf[j - 1]) return(s[j])
  s[j - 1] + (target - yf[j - 1]) * (s[j] - s[j - 1]) / (yf[j] - yf[j - 1])
}

agg <- function(df, by, ...) {
  f <- lapply(by, function(b) df[[b]])
  names(f) <- by
  a <- aggregate(list(...), by = f, FUN = function(v) mean(v, na.rm = TRUE))
  a
}


# ==============================================================================
# 1. SIMULATION INTEGRITY
# ==============================================================================

cat("\n", line, "\n05a DETECTOR REPORT\n", line, "\n", sep = "")

cat("\n1. SIMULATION INTEGRITY\n", subline, "\n", sep = "")

n_eval_files <- length(list.files(checkpoint_dir, pattern = "^eval_cell_"))
n_cal_files  <- length(list.files(checkpoint_dir, pattern = "^cal_cell_"))

pr(data.frame(
  metric = c("Evaluation checkpoints", "Expected", "Calibration checkpoints",
             "Expected", "Summary rows", "Environments", "Seed groups",
             "Min finite MIS rate", "Cells below 99% finite"),
  value = c(n_eval_files, nrow(eval_grid), n_cal_files, nrow(grids$cal_grid),
            nrow(S), length(unique(S$env_id)), length(unique(S$seed_group_id)),
            round(min(S$finite_MIS, na.rm = TRUE), 5),
            sum(S$finite_MIS < 0.99, na.rm = TRUE))))


# ==============================================================================
# 2. COMMON RANDOM NUMBERS
#
# Under linear endogeneity the standardized statistic is provably invariant to
# severity: the confounder is linear in x, so the fitted residual vector is a
# scalar multiple of its severity-0 counterpart and the ratio cancels. If CRN
# holds, mean_MIS is therefore constant across severity within a seed group.
# This is a stronger check than comparing seeds directly.
# ==============================================================================

cat("\n2. COMMON RANDOM NUMBER INTEGRITY (via endogeneity invariance)\n",
    subline, "\n", sep = "")

endo <- P[P$scenario == "endogeneity", ]

crn <- do.call(rbind, lapply(split(endo, endo$seed_group_id), function(g)
  data.frame(seed_group_id = g$seed_group_id[1],
             levels = nrow(g),
             max_abs_dev = max(abs(g$mean_MIS - g$mean_MIS[1]), na.rm = TRUE),
             stringsAsFactors = FALSE)))

pr(data.frame(
  matched_groups = nrow(crn),
  groups_with_10_levels = sum(crn$levels == 10L),
  max_deviation_across_severity = max(crn$max_abs_dev, na.rm = TRUE),
  CRN_pass = max(crn$max_abs_dev, na.rm = TRUE) < 1e-8))

cat("\nA nonzero deviation means either CRN broke or the invariance argument\n")
cat("does not hold as stated. Both are worth knowing before drafting.\n")


# ==============================================================================
# 3. EMPIRICAL SIZE
#
# The number frozen 05 could not produce. Cutoffs come from a disjoint seed
# stream, so this is a real false-positive rate rather than an identity.
# ==============================================================================

cat("\n3. EMPIRICAL SIZE AT SEVERITY 0\n", subline, "\n", sep = "")

sz <- W[W$severity_index == 1L, ]

cat("Pooled over all environments and scenarios, by k:\n\n")
pr(agg(sz, c("k_fraction"),
       size_MIS = sz$rej_MIS, size_Cook = sz$rej_Cook,
       size_DFBETAS = sz$rej_DFBETAS, size_Leverage = sz$rej_Leverage,
       size_RESET_cal = sz$rej_RESET_cal, size_BP_cal = sz$rej_BP_cal,
       size_RESET_nom = sz$rej_RESET_nom, size_BP_nom = sz$rej_BP_nom))

szp <- sz[sz$k_fraction == PRIMARY_K, ]

cat("\nBy error distribution, k = 5%:\n\n")
pr(agg(szp, c("error_type"),
       size_MIS = szp$rej_MIS, size_Cook = szp$rej_Cook,
       size_RESET_cal = szp$rej_RESET_cal, size_BP_cal = szp$rej_BP_cal))

cat("\nBy X distribution and n, k = 5%:\n\n")
pr(agg(szp, c("x_type", "n"), size_MIS = szp$rej_MIS))

cat("\nDispersion of MIS size across the ", nrow(szp), " severity-0 cells:\n",
    sep = "")
pr(data.frame(
  mean = mean(szp$rej_MIS, na.rm = TRUE),
  sd = stats::sd(szp$rej_MIS, na.rm = TRUE),
  p05 = stats::quantile(szp$rej_MIS, 0.05, na.rm = TRUE, names = FALSE),
  p95 = stats::quantile(szp$rej_MIS, 0.95, na.rm = TRUE, names = FALSE),
  exactly_0.05 = sum(abs(szp$rej_MIS - 0.05) < 1e-12)))

cat("\nAny cell exactly 0.0500 would indicate residual circularity.\n")
cat("Nominal vs calibrated columns show how much the size correction was\n")
cat("doing for RESET and BP.\n")


# ==============================================================================
# 4. NULL SCALING IN n
# ==============================================================================

cat("\n4. NULL LEVEL SCALES AS sqrt(n)\n", subline, "\n", sep = "")

nl <- sz[sz$k_fraction == PRIMARY_K, ]
nl$mean_over_sqrt_n <- nl$mean_MIS / sqrt(nl$n)

pr(agg(nl, c("x_type", "n"),
       null_mean_MIS = nl$mean_MIS, mean_over_sqrt_n = nl$mean_over_sqrt_n))

cat("\nlog-log slope of null mean MIS on n (0.5 => sqrt(n)):\n")
for (xt in unique(nl$x_type)) {
  g <- agg(nl[nl$x_type == xt, ], c("n"), m = nl$mean_MIS[nl$x_type == xt])
  if (nrow(g) > 1L) {
    sl <- stats::coef(stats::lm(log(m) ~ log(n), data = g))[2]
    cat("  ", xt, ": ", round(unname(sl), 4), "\n", sep = "")
  }
}

cat("\nThis is why a fixed critical value is meaningless and calibration must\n")
cat("be per-environment.\n")


# ==============================================================================
# 5. DETECTION POWER BY SCENARIO
# ==============================================================================

cat("\n5. DETECTION POWER BY SCENARIO -- k = 5%\n", subline, "\n", sep = "")

pw <- P[P$severity_index > 1L, ]

cat("Mean power over all positive severities:\n\n")
pr(agg(pw, c("scenario"),
       MIS = pw$rej_MIS, Cook = pw$rej_Cook, DFBETAS = pw$rej_DFBETAS,
       Leverage = pw$rej_Leverage, RESET = pw$rej_RESET_cal,
       BP = pw$rej_BP_cal))

mx <- P[P$severity_index == 10L, ]

cat("\nPower at maximum severity:\n\n")
pr(agg(mx, c("scenario"),
       MIS = mx$rej_MIS, Cook = mx$rej_Cook, DFBETAS = mx$rej_DFBETAS,
       Leverage = mx$rej_Leverage, RESET = mx$rej_RESET_cal,
       BP = mx$rej_BP_cal))


# ==============================================================================
# 6. SEVERITY RESPONSE CURVE
# ==============================================================================

cat("\n6. MIS SEVERITY RESPONSE -- k = 5%\n", subline, "\n", sep = "")

pr(agg(P, c("scenario", "severity_index"),
       severity = P$severity_target, mean_MIS = P$mean_MIS,
       MIS_power = P$rej_MIS))

cat("\nWatch for saturation: where mean_MIS stops rising, monotonicity claims\n")
cat("rest on movement in the third decimal even though power is at 1.\n")


# ==============================================================================
# 7. MIS VERSUS RIVALS -- PAIRED NET ADVANTAGE
#
# Net advantage = P(MIS rejects, rival does not) - P(rival rejects, MIS does
# not), computed within draw. Positive favours MIS.
# ==============================================================================

cat("\n7. PAIRED NET ADVANTAGE OF MIS -- k = 5%, positive severities\n",
    subline, "\n", sep = "")

pr(agg(pw, c("scenario"),
       vs_Cook = pw$net_vs_Cook, vs_DFBETAS = pw$net_vs_DFBETAS,
       vs_Leverage = pw$net_vs_Leverage, vs_RESET = pw$net_vs_RESET,
       vs_BP = pw$net_vs_BP))

cat("\nExpect MIS to dominate the influence diagnostics and to lose to each\n")
cat("specification test in its own specialty. The defensible claim is that\n")
cat("MIS is never blind, not that it always wins.\n")


# ==============================================================================
# 8. 80% DETECTION BOUNDARY
# ==============================================================================

cat("\n8. ISOTONIC 80% DETECTION BOUNDARY -- k = 5%\n", subline, "\n", sep = "")

bd_rows <- list()

for (kf in k_grid) {
  
  Pk <- W[W$k_fraction == kf, ]
  
  for (key in unique(paste(Pk$env_id, Pk$scenario, sep = "|"))) {
    parts <- strsplit(key, "|", fixed = TRUE)[[1]]
    g <- Pk[Pk$env_id == as.integer(parts[1]) & Pk$scenario == parts[2], ]
    g <- g[order(g$severity_target), ]
    if (nrow(g) < 3L) next
    bd_rows[[length(bd_rows) + 1L]] <- data.frame(
      env_id = g$env_id[1], n = g$n[1], x_type = g$x_type[1],
      error_type = g$error_type[1], scenario = g$scenario[1],
      k_fraction = kf,
      b_MIS = boundary80(g$severity_target, g$rej_MIS),
      b_Cook = boundary80(g$severity_target, g$rej_Cook),
      b_DFBETAS = boundary80(g$severity_target, g$rej_DFBETAS),
      b_Leverage = boundary80(g$severity_target, g$rej_Leverage),
      b_RESET = boundary80(g$severity_target, g$rej_RESET_cal),
      b_BP = boundary80(g$severity_target, g$rej_BP_cal),
      stringsAsFactors = FALSE)
  }
}

bd_all <- do.call(rbind, bd_rows)
safe_save_rds(bd_all, bound_path)

bd <- bd_all[bd_all$k_fraction == PRIMARY_K, ]

reach <- do.call(rbind, lapply(split(bd, bd$scenario), function(g)
  data.frame(scenario = g$scenario[1], environments = nrow(g),
             reach_MIS = mean(is.finite(g$b_MIS)),
             reach_Cook = mean(is.finite(g$b_Cook)),
             reach_DFBETAS = mean(is.finite(g$b_DFBETAS)),
             reach_Leverage = mean(is.finite(g$b_Leverage)),
             reach_RESET = mean(is.finite(g$b_RESET)),
             reach_BP = mean(is.finite(g$b_BP)),
             stringsAsFactors = FALSE)))

cat("Reach rate -- proportion of environments attaining 80% power:\n\n")
pr(reach)

cat("\nMedian boundary among environments that reached 80%\n")
cat("(conditional on reaching; read alongside the reach rate above):\n\n")
pr(do.call(rbind, lapply(split(bd, bd$scenario), function(g)
  data.frame(scenario = g$scenario[1],
             MIS = stats::median(g$b_MIS, na.rm = TRUE),
             n_MIS = sum(is.finite(g$b_MIS)),
             RESET = stats::median(g$b_RESET, na.rm = TRUE),
             n_RESET = sum(is.finite(g$b_RESET)),
             BP = stats::median(g$b_BP, na.rm = TRUE),
             n_BP = sum(is.finite(g$b_BP)),
             stringsAsFactors = FALSE))))

cat("\nMIS boundary by n (median over environments reaching 80%):\n\n")
pr(do.call(rbind, lapply(split(bd, list(bd$scenario, bd$n), drop = TRUE),
                         function(g) data.frame(scenario = g$scenario[1], n = g$n[1],
                                                median_boundary = stats::median(g$b_MIS, na.rm = TRUE),
                                                reach_rate = mean(is.finite(g$b_MIS)),
                                                stringsAsFactors = FALSE))))


# ==============================================================================
# 9. EFFECT OF THE DELETION BUDGET k
#
# Frozen 05 showed severity correlation collapsing between k = 5% and 10%.
# With a valid calibration this can finally be read as a property of the
# estimator rather than an artifact of the cutoff.
# ==============================================================================

cat("\n9. EFFECT OF DELETION BUDGET k\n", subline, "\n", sep = "")

kk <- W[W$severity_index > 1L, ]

pr(agg(kk, c("scenario", "k_fraction"), MIS_power = kk$rej_MIS))

cat("\nSize by k (severity 0):\n\n")
pr(agg(sz, c("k_fraction"), size_MIS = sz$rej_MIS))


# ==============================================================================
# 10. RESPECIFICATION -- SPECIFICITY OF THE DETECTOR
# ==============================================================================

cat("\n10. WRONG VERSUS CORRECT SPECIFICATION -- k = 5%\n", subline, "\n",
    sep = "")

cs <- S[S$k_fraction == PRIMARY_K &
          S$scenario %in% unique(S$scenario[S$model_state == "correct"]), ]

wide <- merge(
  cs[cs$model_state == "wrong",
     c("cell_id", "scenario", "n", "x_type", "error_type",
       "severity_index", "severity_target", "mean_MIS", "rej_MIS")],
  cs[cs$model_state == "correct",
     c("cell_id", "mean_MIS", "rej_MIS")],
  by = "cell_id", suffixes = c("_wrong", "_correct"))

wide$reduction <- wide$mean_MIS_wrong - wide$mean_MIS_correct

pos <- wide[wide$severity_index > 1L, ]

pr(agg(pos, c("scenario"),
       mean_wrong = pos$mean_MIS_wrong, mean_correct = pos$mean_MIS_correct,
       mean_reduction = pos$reduction,
       prop_correct_lower = as.numeric(pos$reduction > 0),
       power_wrong = pos$rej_MIS_wrong, power_correct = pos$rej_MIS_correct))

cat("\npower_correct is the false-positive rate after the model is fixed and\n")
cat("should sit near the nominal size. mean_reduction > 0 means correcting\n")
cat("the specification removed the MIS signal.\n")

cat("\nBy X distribution -- pilot 1 found the threshold anomaly was confined\n")
cat("to contaminated X (over-specification inflation under heavy tails):\n\n")
pr(agg(pos, c("scenario", "x_type"), mean_reduction = pos$reduction))


# ==============================================================================
# 11. LOCALISATION
# ==============================================================================

cat("\n11. SUBGROUP RECOVERY -- k = 5%\n", subline, "\n", sep = "")

loc <- P[P$severity_index > 1L & is.finite(P$mean_lift), ]

if (nrow(loc) > 0L) {
  pr(agg(loc, c("scenario"),
         affected_fraction = loc$mean_affected,
         precision = loc$mean_precision, recall = loc$mean_recall,
         lift = loc$mean_lift))
  cat("\nRecall is bounded above by k / affected_fraction, so precision and\n")
  cat("lift are the interpretable columns. Neither RESET nor BP produces\n")
  cat("any of these -- localisation is what MIS adds.\n")
} else {
  cat("No localisation rows (no scenario recorded affected_idx).\n")
}


# ==============================================================================
# 12. ENDOGENEITY INVARIANCE AND ITS SCOPE
# ==============================================================================

cat("\n12. ENDOGENEITY: LINEAR VERSUS NONLINEAR CONFOUNDING -- k = 5%\n",
    subline, "\n", sep = "")

en <- P[P$scenario %in% c("endogeneity", "endogeneity_nl"), ]

pr(agg(en, c("scenario", "severity_index"),
       severity = en$severity_target, mean_MIS = en$mean_MIS,
       MIS_power = en$rej_MIS, Cook_power = en$rej_Cook,
       RESET_power = en$rej_RESET_cal))

cat("\nLinear endogeneity should be exactly flat for every influence-based\n")
cat("diagnostic. If endogeneity_nl rises, the invariance result is specific\n")
cat("to linear-in-x confounding and must be stated that way.\n")


# ==============================================================================
# 13. ENVIRONMENT HETEROGENEITY
# ==============================================================================

cat("\n13. ROBUSTNESS ACROSS ENVIRONMENTS -- k = 5%, positive severities\n",
    subline, "\n", sep = "")

cat("By error distribution:\n\n")
pr(agg(pw, c("error_type"),
       MIS = pw$rej_MIS, Cook = pw$rej_Cook, RESET = pw$rej_RESET_cal,
       BP = pw$rej_BP_cal))

cat("\nBy X distribution:\n\n")
pr(agg(pw, c("x_type"),
       MIS = pw$rej_MIS, Cook = pw$rej_Cook, RESET = pw$rej_RESET_cal,
       BP = pw$rej_BP_cal))

cat("\nBy n:\n\n")
pr(agg(pw, c("n"),
       MIS = pw$rej_MIS, Cook = pw$rej_Cook, RESET = pw$rej_RESET_cal,
       BP = pw$rej_BP_cal))

cat("\nWorst-case MIS power by scenario across the 96 environments:\n\n")
pr(do.call(rbind, lapply(split(mx, mx$scenario), function(g)
  data.frame(scenario = g$scenario[1],
             min_power = min(g$rej_MIS, na.rm = TRUE),
             p10 = stats::quantile(g$rej_MIS, 0.10, na.rm = TRUE, names = FALSE),
             median = stats::median(g$rej_MIS, na.rm = TRUE),
             max_power = max(g$rej_MIS, na.rm = TRUE),
             stringsAsFactors = FALSE))))


# ==============================================================================
# VERDICT
# ==============================================================================

cat("\n", line, "\nAUTOMATED VERDICT\n", line, "\n", sep = "")

size_ok <- abs(mean(szp$rej_MIS, na.rm = TRUE) - 0.05) < 0.015
degen   <- sum(abs(szp$rej_MIS - 0.05) < 1e-12)

cat("Pooled MIS size at k=5%: ", round(mean(szp$rej_MIS, na.rm = TRUE), 4),
    "  -> ", if (size_ok) "PASS" else "CHECK", "\n", sep = "")
cat("Cells exactly 0.0500 (circularity flag): ", degen, "\n", sep = "")
cat("CRN via endogeneity invariance: ",
    if (max(crn$max_abs_dev, na.rm = TRUE) < 1e-8) "PASS" else "CHECK",
    "\n", sep = "")
cat("Missing evaluation checkpoints: ",
    nrow(eval_grid) - n_eval_files, "\n", sep = "")
cat("Cells below 99% finite MIS: ", sum(S$finite_MIS < 0.99, na.rm = TRUE),
    "\n", sep = "")

cat("\nSaved: ", summary_path, "\n       ", bound_path, "\n", sep = "")
cat(line, "\n", sep = "")