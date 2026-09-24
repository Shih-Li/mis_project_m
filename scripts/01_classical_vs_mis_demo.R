# ==============================================================================
# Script: /script/01_classical_vs_mis_demo.R
# Purpose: Demonstrate classical LOO failure vs MIS under masking.
# Inputs: None (Synthetic Data Generation Process defined internally)
# Outputs: ../output/fig1_classical_vs_mis_demo.pdf&png
# ==============================================================================

# 1. Setup
script_file <- normalizePath(
  sys.frame(1)$ofile,
  winslash = "/",
  mustWork = TRUE
)

script_dir <- dirname(script_file)

project_dir <- normalizePath(
  file.path(script_dir, ".."),
  winslash = "/",
  mustWork = TRUE
)

output_dir <- file.path(
  project_dir,
  "output"
)

source(
  file.path(
    project_dir,
    "R",
    "diagnostics_classical.R"
  )
)

source(
  file.path(
    project_dir,
    "R",
    "dinkelbach_topk.R"
  )
)

cat(
  "\nScript directory:",
  script_dir,
  "\n"
)

cat(
  "Project directory:",
  project_dir,
  "\n"
)

cat(
  "Output directory:",
  output_dir,
  "\n"
)

cat(
  "Output directory exists:",
  dir.exists(output_dir),
  "\n\n"
)

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

set.seed(123)

# 2. DGP
# ------------------------------------------------------------------
# Total N = 400
# Bad-leverage set = observations 1:10
# 10 / 400 = 2.5%
#
# Three highlighted contamination types:
#   1:10  = bad leverage
#   11:20 = good leverage
#   21:30 = response outliers
#   31:400 = regular observations
# ------------------------------------------------------------------

n <- 400L

# cluster size
k_cluster <- 10L

# total bad-leverage coalition size = 20
k_bad <- 2L * k_cluster

# Common deletion budget for all methods
k_select <- k_bad

# MIS uses the same deletion budget
k_mis <- k_select

stopifnot(k_bad / n == 0.05)

ids_bad1      <- 1:10
ids_bad2      <- 11:20
ids_bad       <- c(ids_bad1, ids_bad2)

ids_good      <- 21:30
ids_response  <- 31:40
ids_regular   <- 41:n

# ------------------------------------------------------------
# Bad leverage: observations 1:10
#
# Moderately displaced X values, all pushing the slope in the
# same direction. Individually they are kept below the usual
# classical thresholds, but jointly they form the influential set.
# ------------------------------------------------------------
# Bad 1: upper bad-leverage cluster
x_bad1 <- rnorm(k_cluster, mean = 2.35, sd = 0.05)
y_bad1 <- rnorm(k_cluster, mean = -0.15, sd = 0.05)

# Bad 2: lower bad-leverage cluster
x_bad2 <- rnorm(k_cluster, mean = 2.75, sd = 0.05)
y_bad2 <- rnorm(k_cluster, mean = -1.00, sd = 0.05)

# ------------------------------------------------------------
# Good leverage: observations 11:20
#
# Very large |X| but exactly on the true regression line.
# These should be obvious to leverage but should not distort beta.
# ------------------------------------------------------------
x_good <- rnorm(k_cluster, mean = 8.8, sd = 0.10)
y_good <- 1.5 * x_good + rnorm(k_cluster, sd = 0.05)

# ------------------------------------------------------------
# Response outliers: observations 21:30
#
# X remains near the centre of the design, while Y is moved
# strongly away from the regression line.
# These are designed to be obvious to Cook's distance.
# ------------------------------------------------------------
x_response <- rnorm(k_cluster, mean = 0.0, sd = 0.08)
y_response <- rnorm(k_cluster, mean = 20.0, sd = 0.20)

# ------------------------------------------------------------
# Regular observations: 31:400
# ------------------------------------------------------------
n_regular <- n - 4L * k_cluster

x_regular <- rnorm(n_regular)
y_regular <- 1.5 * x_regular +
  rnorm(n_regular, sd = 0.5)

df <- data.frame(
  id = seq_len(n),
  
  x = c(
    x_bad1,
    x_bad2,
    x_good,
    x_response,
    x_regular
  ),
  
  y = c(
    y_bad1,
    y_bad2,
    y_good,
    y_response,
    y_regular
  ),
  
  truth = factor(
    c(
      rep("Bad leverage 1",             k_cluster),
      rep("Bad leverage 2",             k_cluster),
      rep("Good leverage",     k_cluster),
      rep("Response outlier",  k_cluster),
      rep("Regular", n_regular)
    ),
    levels = c(
      "Regular",
      "Good leverage",
      "Bad leverage 1",
      "Bad leverage 2",
      "Response outlier"
    )
  )
)

mdl_base <- lm(y ~ x, data = df)

# ------------------------------------------------------------------
# Quantify the effect of the planted bad-leverage coalition
# ------------------------------------------------------------------

mdl_without_bad <- lm(
  y ~ x,
  data = df[-ids_bad, ]
)

beta_full <- unname(
  coef(mdl_base)["x"]
)

beta_without_bad <- unname(
  coef(mdl_without_bad)["x"]
)

joint_beta_shift <- abs(
  beta_without_bad - beta_full
)

joint_beta_shift_pct <- 100 *
  joint_beta_shift / 1.5

# ------------------------------------------------------------------
# Measure slope error relative to the true coefficient
# ------------------------------------------------------------------

beta_true <- 1.5

# Error in the contaminated OLS estimate
bias_before_mis <- abs(
  beta_full - beta_true
)

# Error after removing the MIS / true bad-leverage coalition
bias_after_mis <- abs(
  beta_without_bad - beta_true
)

# Percentage reduction in slope error after MIS removal
bias_reduction_pct <- 100 * (
  1 - bias_after_mis / bias_before_mis
)

cat(
  sprintf(
    paste0(
      "Slope error before MIS removal: %.3f\n",
      "Slope error after MIS removal: %.3f\n",
      "Slope error reduction: %.1f%%\n\n"
    ),
    bias_before_mis,
    bias_after_mis,
    bias_reduction_pct
  )
)


# ------------------------------------------------------------------
# Leave-one-out effects for the 10 planted bad-leverage observations
#
# These are the effects classical one-at-a-time diagnostics are
# trying to detect. Under masking, each one looks small while the
# whole coalition has a much larger joint effect.
# ------------------------------------------------------------------

beta_loo_bad <- vapply(
  ids_bad,
  function(i) {
    unname(
      coef(
        lm(
          y ~ x,
          data = df[-i, ]
        )
      )["x"]
    )
  },
  numeric(1)
)

loo_beta_shift <- abs(
  beta_loo_bad - beta_full
)

max_loo_shift <- max(
  loo_beta_shift
)

masking_ratio <- joint_beta_shift /
  max_loo_shift

cat(
  sprintf(
    paste0(
      "Slope with all observations: %.3f\n",
      "Slope after removing bad-leverage set: %.3f\n",
      "Joint slope change: %.3f\n",
      "Joint change as %% of true slope: %.1f%%\n",
      "Largest single-deletion change: %.4f\n",
      "Joint / largest single effect: %.1fx\n\n"
    ),
    beta_full,
    beta_without_bad,
    joint_beta_shift,
    joint_beta_shift_pct,
    max_loo_shift,
    masking_ratio
  )
)

# 3. Detection
# ------------------------------------------------------------------
# Fixed-budget comparison:
#
# Every method is allowed to remove exactly k_select = 20 observations.
# Classical diagnostics therefore contribute their top-20 ranked points,
# while MIS jointly selects its most influential 20-point subset.
#
# This makes the comparison conditional on the same deletion budget.
# ------------------------------------------------------------------

ids_cooks <- get_classical_set(
  mdl_base,
  target_var = "x",
  k = k_select,
  metric = "cooks_d"
)

ids_lev <- get_classical_set(
  mdl_base,
  target_var = "x",
  k = k_select,
  metric = "leverage"
)

ids_dfb <- get_classical_set(
  mdl_base,
  target_var = "x",
  k = k_select,
  metric = "dfbetas_target"
)

# Exact fixed-k MIS.
# sign = -1 targets the direction generated by the planted
# bad-leverage coalition.
ids_mis <- dinkelbach_topk_lm(
  mdl_base,
  pos  = 2L,
  sign = -1L,
  k    = k_select
)

flagged <- list(
  cooks = ids_cooks,
  lev   = ids_lev,
  dfb   = ids_dfb,
  mis   = ids_mis
)

# ------------------------------------------------------------------
# Demo validation
# ------------------------------------------------------------------

cooks_bad_hits <- intersect(
  ids_cooks,
  ids_bad
)

lev_bad_hits <- intersect(
  ids_lev,
  ids_bad
)

dfb_bad_hits <- intersect(
  ids_dfb,
  ids_bad
)

mis_bad_hits <- intersect(
  ids_mis,
  ids_bad
)

slope_after_delete <- function(idx) {
  
  if (length(idx) == 0L) {
    return(beta_full)
  }
  
  unname(
    coef(
      lm(
        y ~ x,
        data = df[-idx, , drop = FALSE]
      )
    )["x"]
  )
}

cat(
  sprintf(
    "\nBad-leverage contamination: %d / %d = %.1f%%\n",
    k_bad,
    n,
    100 * k_bad / n
  )
)

cat(
  sprintf(
    paste0(
      "Bad-leverage recovery under common k = %d:\n",
      "  Cook's D : %d / %d\n",
      "  Leverage : %d / %d\n",
      "  DFBETAS  : %d / %d\n",
      "  MIS      : %d / %d\n\n"
    ),
    k_select,
    length(cooks_bad_hits), k_bad,
    length(lev_bad_hits),   k_bad,
    length(dfb_bad_hits),   k_bad,
    length(mis_bad_hits),   k_bad
  )
)

# Report whether the intended demonstration was obtained.
demo_ok <- (
  length(ids_cooks) == k_select &&
    length(ids_lev) == k_select &&
    length(ids_dfb) == k_select &&
    length(ids_mis) == k_select
)

stopifnot(
  length(ids_cooks) == k_select,
  length(ids_lev)   == k_select,
  length(ids_dfb)   == k_select,
  length(ids_mis)   == k_select
)

cat(
  "Demo validation:",
  if (demo_ok) "PASS" else "WARNING - target pattern not fully obtained",
  "\n\n"
)

# ------------------------------------------------------------------
# Coefficient comparison for Panel 2
# ------------------------------------------------------------------

coef_compare <- data.frame(
  set = c(
    "Contaminated\nOLS",
    "Cook's D",
    "Leverage",
    "DFBETAS",
    "MIS",
    "True bad\nset"
  ),
  
  beta = c(
    beta_full,
    slope_after_delete(ids_cooks),
    slope_after_delete(ids_lev),
    slope_after_delete(ids_dfb),
    slope_after_delete(ids_mis),
    slope_after_delete(ids_bad)
  )
)

coef_compare$error <- abs(
  coef_compare$beta - beta_true
)

recovery_compare <- c(
  NA,
  length(cooks_bad_hits),
  length(lev_bad_hits),
  length(dfb_bad_hits),
  length(mis_bad_hits),
  k_bad
)

# Expected pedagogical structure:
#
#   Good leverage      -> classical leverage can identify it.
#   Response outlier   -> Cook's distance can identify it.
#   Bad leverage 1:10  -> classical LOO diagnostics miss the coalition.
#                         Exact fixed-k MIS recovers 1:10 jointly.
#
# The important change relative to the current script is that the
# three colours describe the TRUE contamination topology, not the
# diagnostic method. Detection itself is shown with black outlines.
#
# Also, classical diagnostics now use their normal thresholds rather
# than forcing k = 10. MIS alone is given the known deletion budget
# k = 10 = 2.5% of N.

# 4. Generate Output
# ------------------------------------------------------------------
# Truth colours
#
# There are three substantive colours:
#   blue       = good leverage
#   red/orange = bad leverage
#   purple     = response outlier
#
# Regular observations stay neutral grey.
# ------------------------------------------------------------------

truth_cols <- c(
  "Regular"          = "gray80",
  "Good leverage"    = "#0072B2",
  "Bad leverage 1"            = "#D55E00",
  "Bad leverage 2"            = "#E69F00",
  "Response outlier" = "#CC79A7"
)

point_cols <- truth_cols[as.character(df$truth)]

draw_base <- function(main_title, xlim = NULL, ylim = NULL) {
  
  if (is.null(xlim)) {
    xlim <- range(df$x) + c(-0.5, 0.5)
  }
  
  if (is.null(ylim)) {
    ylim <- range(df$y) + c(-1, 1)
  }
  
  plot(
    df$x,
    df$y,
    main = main_title,
    xlab = "X",
    ylab = "Y",
    pch = 16,
    col = point_cols,
    cex = 0.9,
    xlim = xlim,
    ylim = ylim
  )
  
  # True DGP line
  abline(
    a = 0,
    b = 1.5,
    lty = 2,
    lwd = 1.5,
    col = "gray35"
  )
}


draw_panels <- function() {
  
  par(
    mfrow = c(1, 2),
    mar = c(5, 5, 4, 1.5)
  )
  
  
  # ============================================================
  # Panel 1
  # Cluster geometry
  # ============================================================
  
  draw_base(
    "Contamination Topology (Each Cluster Has n = 10)"
  )
  
  # Highlight the contamination clusters
  points(
    df$x[ids_good],
    df$y[ids_good],
    pch = 21,
    bg = truth_cols["Good leverage"],
    col = "black",
    cex = 1.5,
    lwd = 1.1
  )
  
  points(
    df$x[ids_bad1],
    df$y[ids_bad1],
    pch = 21,
    bg = truth_cols["Bad leverage 1"],
    col = "black",
    cex = 1.5,
    lwd = 1.1
  )
  
  points(
    df$x[ids_bad2],
    df$y[ids_bad2],
    pch = 21,
    bg = truth_cols["Bad leverage 2"],
    col = "black",
    cex = 1.5,
    lwd = 1.1
  )
  
  points(
    df$x[ids_response],
    df$y[ids_response],
    pch = 21,
    bg = truth_cols["Response outlier"],
    col = "black",
    cex = 1.5,
    lwd = 1.1
  )
  
  legend(
    "topleft",
    legend = c(
      "Regular",
      "Good leverage",
      "Bad leverage 1",
      "Bad leverage 2",
      "Response outlier"
    ),
    pch = 16,
    col = truth_cols[
      c(
        "Regular",
        "Good leverage",
        "Bad leverage 1",
        "Bad leverage 2",
        "Response outlier"
      )
    ],
    pt.cex = 1.1,
    cex = 0.82,
    bg = "white"
  )
  
  # ============================================================
  # Panel 2
  # Coefficient change after deleting different 10-point sets
  # ============================================================
  
  bp <- barplot(
    coef_compare$error,
    names.arg = coef_compare$set,
    las = 1,
    ylim = c(
      0,
      max(coef_compare$error) * 1.25
    ),
    col = c(
      "gray60",
      "gray75",
      "gray75",
      "gray75",
      truth_cols["Bad leverage"],
      "gray90"
    ),
    border = "gray30",
    ylab = "Absolute slope error",
    main = "Same 20-Point Budget: Which Method Finds the Harmful Points?",
    cex.names = 0.78
  )
  
  bar_labels <- ifelse(
    is.na(recovery_compare),
    sprintf(
      "error = %.3f",
      coef_compare$error
    ),
    sprintf(
      "error = %.3f\nbad = %d/%d",
      coef_compare$error,
      recovery_compare,
      k_bad
    )
  )
  
  text(
    x = bp,
    y = coef_compare$error,
    labels = bar_labels,
    pos = 3,
    cex = 0.75
  )
  
  text(
    x = bp,
    y = coef_compare$beta,
    labels = sprintf("%.3f", coef_compare$beta),
    pos = 3,
    cex = 0.82
  )
  
  text(
    x = max(bp),
    y = 1.5,
    labels = "True slope = 1.500",
    pos = 3,
    cex = 0.82,
    col = "gray35"
  )
  
  mtext(
    sprintf(
      "Every method removes %d observations; bars show remaining slope error",
      k_select
    ),
    side = 3,
    line = 0.4,
    cex = 0.80
  )
}

## ------------------------------------------------------------------
# Export
# ------------------------------------------------------------------

pdf_file <- file.path(
  output_dir,
  "fig1_classical_vs_mis_demo.pdf"
)

png_file <- file.path(
  output_dir,
  "fig1_classical_vs_mis_demo.png"
)

cat(
  "Writing PDF to:\n",
  pdf_file,
  "\n\n"
)

pdf(
  pdf_file,
  width = 14,
  height = 6
)

draw_panels()

dev.off()


cat(
  "Writing PNG to:\n",
  png_file,
  "\n\n"
)

png(
  png_file,
  width = 14,
  height = 6,
  units = "in",
  res = 300
)

draw_panels()

dev.off()


cat(
  "PDF exists:",
  file.exists(pdf_file),
  "\n"
)

cat(
  "PNG exists:",
  file.exists(png_file),
  "\n\n"
)

cat(
  sprintf(
    paste0(
      "Simulation complete.\n",
      "N = %d; true bad-leverage set = %d observations (%.1f%%).\n",
      "Common deletion budget = %d.\n",
      "Cook's D recovery = %d/%d.\n",
      "Leverage recovery = %d/%d.\n",
      "DFBETAS recovery = %d/%d.\n",
      "MIS recovery = %d/%d.\n"
    ),
    n,
    k_bad,
    100 * k_bad / n,
    k_select,
    length(cooks_bad_hits), k_bad,
    length(lev_bad_hits),   k_bad,
    length(dfb_bad_hits),   k_bad,
    length(mis_bad_hits),   k_bad
  )
)