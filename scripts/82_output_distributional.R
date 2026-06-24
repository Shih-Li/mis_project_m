# ==============================================================================
# File: scripts/82_output_distributional.R
# Purpose:
#   Generate publication-ready figures and tables for:
#     1. Script 02: broad distributional robustness experiment
#     2. Script 02b: good-leverage mechanism diagnostic
#
# Inputs:
#   output/02_final_distributions.rds
#   output/02b_good_leverage_diagnostic.rds
#
# Output structure:
#   output/02_distributional/
#     figures/main/          Main-paper PDF figures
#     figures/supplement/    Supplementary PDF figures
#     tables/main/           Complete LaTeX table environments + CSV companions
#     tables/supplement/     Supplementary LaTeX tables + CSV companions
#     data/                  Data used to construct figures and tables
#     diagnostics/           Session information and output manifest
#
# LaTeX requirements for generated tables:
#   \usepackage{float}
#   \usepackage{graphicx}
#   \usepackage{booktabs}
#
# Design rules:
#   - No title or subtitle inside plots; captions are handled in LaTeX.
#   - Blue-orange diverging heatmaps: orange = MIS worse, blue = MIS better.
#   - Sequential quantities use blue, never red-green.
#   - Main tables are complete \begin{table}[H] environments.
#   - Compact tables use 0.7\columnwidth; wider tables use a larger width only
#     when needed for readable labels and indicators.
# ==============================================================================

# ==============================================================================
# 0. Packages, paths, and configuration
# ==============================================================================

required_packages <- c("dplyr", "tidyr", "ggplot2", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running scripts/82_output_distributional.R."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

resolve_project_root <- function() {
  candidates <- unique(c(
    normalizePath(".",  winslash = "/", mustWork = FALSE),
    normalizePath("..", winslash = "/", mustWork = FALSE)
  ))
  valid <- candidates[
    dir.exists(file.path(candidates, "R")) &
      dir.exists(file.path(candidates, "scripts"))
  ]
  if (length(valid) == 0L) {
    stop(
      "Cannot locate project root. Run this script from the repository root ",
      "or from the scripts/ directory."
    )
  }
  valid[[1L]]
}

project_root <- resolve_project_root()
input_02     <- file.path(project_root, "output", "02_final_distributions.rds")
input_02b    <- file.path(project_root, "output", "02b_good_leverage_diagnostic.rds")
output_root  <- file.path(project_root, "output", "02_distributional")

fig_main_dir <- file.path(output_root, "figures", "main")
fig_supp_dir <- file.path(output_root, "figures", "supplement")
tab_main_dir <- file.path(output_root, "tables", "main")
tab_supp_dir <- file.path(output_root, "tables", "supplement")
data_dir     <- file.path(output_root, "data")
diag_dir     <- file.path(output_root, "diagnostics")

for (path in c(
  fig_main_dir, fig_supp_dir, tab_main_dir, tab_supp_dir, data_dir, diag_dir
)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(input_02)) {
  stop("Missing input file: ", input_02)
}
if (!file.exists(input_02b)) {
  stop("Missing input file: ", input_02b)
}

# PDF is the publication format. Set TRUE only when raster previews are useful.
SAVE_PNG_PREVIEWS <- FALSE
PNG_DPI <- 320

# Color-blind-conscious palette.
COL_ORANGE <- "#E69F00"
COL_BLUE   <- "#0072B2"
COL_BLUE_DARK <- "#08519C"
COL_NEUTRAL <- "#F7F7F7"
COL_GREY <- "#4D4D4D"

TABLE_WIDTH_COMPACT <- "0.7\\columnwidth"
TABLE_WIDTH_MEDIUM  <- "0.85\\columnwidth"
TABLE_WIDTH_WIDE    <- "\\columnwidth"

# ==============================================================================
# 1. Shared helper functions
# ==============================================================================

safe_mean <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  unname(stats::quantile(x, probs = probability, na.rm = TRUE, names = FALSE))
}

row_max_na <- function(...) {
  values <- cbind(...)
  apply(values, 1L, function(x) {
    if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  })
}

safe_ratio <- function(numerator, denominator) {
  out <- rep(NA_real_, length(numerator))
  valid <- is.finite(numerator) & is.finite(denominator) &
    abs(denominator) > sqrt(.Machine$double.eps)
  out[valid] <- numerator[valid] / denominator[valid]
  out
}

safe_log2_ratio <- function(numerator, denominator, use_absolute = FALSE) {
  if (use_absolute) {
    numerator <- abs(numerator)
    denominator <- abs(denominator)
  }
  ratio <- safe_ratio(numerator, denominator)
  ratio[!is.finite(ratio) | ratio <= 0] <- NA_real_
  log2(ratio)
}

format_contamination <- function(x) {
  scales::percent(x, accuracy = 0.1, trim = TRUE)
}

# Script 02 was designed around five nominal contamination proportions.
# The saved contam_prop equals set_size / n_obs after integer rounding, so
# nominally identical design levels can appear as slightly different values
# (for example, 0.0048 and 0.005). Mapping back to the nominal grid prevents
# duplicate factor labels such as two separate "0.5%" levels.
NOMINAL_CONTAM_LEVELS <- c(0.005, 0.01, 0.025, 0.05, 0.10)

map_to_nominal_contamination <- function(
    x,
    nominal_levels = NOMINAL_CONTAM_LEVELS
) {
  vapply(
    x,
    function(value) {
      if (!is.finite(value)) return(NA_real_)
      nominal_levels[which.min(abs(nominal_levels - value))]
    },
    numeric(1)
  )
}

fmt_pct <- function(x, digits = 1L) {
  ifelse(
    is.na(x),
    "--",
    paste0(formatC(100 * x, format = "f", digits = digits), "%")
  )
}

fmt_pp <- function(x, digits = 1L, signed = TRUE) {
  ifelse(
    is.na(x),
    "--",
    if (signed) {
      sprintf(paste0("%+.", digits, "f pp"), 100 * x)
    } else {
      sprintf(paste0("%.", digits, "f pp"), 100 * x)
    }
  )
}

fmt_num <- function(x, digits = 2L) {
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits))
}

fmt_interval <- function(center, lower, upper, digits = 2L) {
  ifelse(
    is.na(center),
    "--",
    paste0(
      formatC(center, format = "f", digits = digits),
      " [",
      formatC(lower, format = "f", digits = digits),
      ", ",
      formatC(upper, format = "f", digits = digits),
      "]"
    )
  )
}

fmt_pct_interval <- function(center, lower, upper, digits = 1L) {
  ifelse(
    is.na(center),
    "--",
    paste0(
      formatC(100 * center, format = "f", digits = digits),
      "% [",
      formatC(100 * lower, format = "f", digits = digits),
      ", ",
      formatC(100 * upper, format = "f", digits = digits),
      "]"
    )
  )
}

escape_latex <- function(x) {
  replacements <- c(
    "\\" = "\\textbackslash{}",
    "&" = "\\&",
    "%" = "\\%",
    "$" = "\\$",
    "#" = "\\#",
    "_" = "\\_",
    "{" = "\\{",
    "}" = "\\}",
    "~" = "\\textasciitilde{}",
    "^" = "\\textasciicircum{}"
  )
  
  escape_one <- function(value) {
    if (is.na(value)) return("")
    characters <- strsplit(as.character(value), "", fixed = TRUE)[[1L]]
    mapped <- replacements[characters]
    characters[!is.na(mapped)] <- unname(mapped[!is.na(mapped)])
    paste0(characters, collapse = "")
  }
  
  vapply(x, escape_one, character(1), USE.NAMES = FALSE)
}

write_tex_table <- function(
    data,
    tex_path,
    caption,
    label,
    resize_width = TABLE_WIDTH_COMPACT,
    align = NULL,
    csv_path = sub("\\.tex$", ".csv", tex_path),
    placement = "H"
) {
  if (!is.data.frame(data) || ncol(data) == 0L) {
    stop("write_tex_table() requires a non-empty data.frame.")
  }
  
  utils::write.csv(data, csv_path, row.names = FALSE, na = "")
  
  if (is.null(align)) {
    align <- paste0("l", paste(rep("r", max(0L, ncol(data) - 1L)), collapse = ""))
  }
  if (nchar(align) != ncol(data)) {
    stop("Length of LaTeX alignment string must equal the number of columns.")
  }
  
  escaped_names <- escape_latex(names(data))
  escaped_data <- lapply(data, escape_latex)
  escaped_data <- as.data.frame(escaped_data, stringsAsFactors = FALSE)
  
  body_rows <- apply(escaped_data, 1L, function(row) {
    paste0(paste(row, collapse = " & "), " \\\\")
  })
  
  latex <- c(
    paste0("\\begin{table}[", placement, "]"),
    "\\centering",
    paste0("\\caption{", escape_latex(caption), "}"),
    paste0("\\label{", label, "}"),
    paste0("\\resizebox{", resize_width, "}{!}{%"),
    paste0("\\begin{tabular}{@{}", align, "@{}}"),
    "\\toprule",
    paste0(paste(escaped_names, collapse = " & "), " \\\\"),
    "\\midrule",
    body_rows,
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\end{table}"
  )
  
  writeLines(latex, tex_path, useBytes = TRUE)
  invisible(tex_path)
}

save_plot <- function(plot, filename, width, height) {
  pdf_path <- file.path(dirname(filename), paste0(tools::file_path_sans_ext(basename(filename)), ".pdf"))
  ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    limitsize = FALSE
  )
  
  if (isTRUE(SAVE_PNG_PREVIEWS)) {
    png_path <- file.path(dirname(filename), paste0(tools::file_path_sans_ext(basename(filename)), ".png"))
    ggsave(
      filename = png_path,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = PNG_DPI,
      limitsize = FALSE
    )
  }
  invisible(pdf_path)
}

theme_paper <- function(base_size = 10.5, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = base_size + 0.5),
      axis.text = element_text(size = base_size),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.25),
      legend.position = "bottom",
      legend.box = "horizontal",
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      panel.spacing = grid::unit(1.2, "lines"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

theme_heatmap <- function(base_size = 10.5, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = base_size + 0.5),
      axis.text = element_text(size = base_size),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.25),
      legend.position = "bottom",
      strip.text = element_text(size = base_size, face = "bold"),
      panel.grid = element_blank(),
      panel.spacing = grid::unit(1.0, "lines"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

# ==============================================================================
# 2. Load and validate inputs
# ==============================================================================

sim <- readRDS(input_02)
diag_02b <- readRDS(input_02b)

required_02_columns <- c(
  "iter", "n_obs", "x_type", "error_type", "dist_param",
  "outlier_method", "set_size", "contam_prop", "block_count",
  "detection_success", "detect_cooks", "detect_lev", "detect_dfbetas",
  "overlap_mis", "overlap_cooks", "overlap_lev", "overlap_dfbetas",
  "cover_evd", "cover_cooks", "cover_lev", "cover_dfbetas",
  "compute_time", "shape", "scale", "loc", "p_value", "converged"
)
missing_02_columns <- setdiff(required_02_columns, names(sim))
if (length(missing_02_columns) > 0L) {
  stop(
    "02_final_distributions.rds is missing required column(s): ",
    paste(missing_02_columns, collapse = ", ")
  )
}

required_02b_columns <- c(
  "error_type", "overlap_mis", "dfb_detected", "dfb_injected",
  "mean_lev_detected", "mean_lev_injected",
  "mean_res_detected", "mean_res_injected"
)
missing_02b_columns <- setdiff(required_02b_columns, names(diag_02b))
if (length(missing_02b_columns) > 0L) {
  stop(
    "02b_good_leverage_diagnostic.rds is missing required column(s): ",
    paste(missing_02b_columns, collapse = ", ")
  )
}

if (!"mix_prop" %in% names(sim)) {
  warning(
    "mix_prop is not stored in 02_final_distributions.rds. Results can be ",
    "aggregated across mixture settings, but the mixture proportions cannot ",
    "be shown separately. Consider adding mix_prop to R/sim_engine.R in a ",
    "future simulation run."
  )
}

cat(sprintf("Loaded Script 02 data:  %s rows\n", format(nrow(sim), big.mark = ",")))
cat(sprintf("Loaded Script 02b data: %s rows\n", format(nrow(diag_02b), big.mark = ",")))

# ==============================================================================
# 3. Labels and tidy preparation
# ==============================================================================

outlier_levels <- c("none", "vertical_outlier", "good_leverage", "bad_leverage")
outlier_labels <- c(
  "none" = "No contamination",
  "vertical_outlier" = "Vertical outliers",
  "good_leverage" = "Good leverage",
  "bad_leverage" = "Bad leverage"
)

error_levels <- c(
  "normal", "beta_logistic", "mixed_normal", "skewed_t",
  "contaminated", "golm", "pareto", "gpd"
)
error_labels <- c(
  "normal" = "Normal",
  "beta_logistic" = "Beta-logistic",
  "mixed_normal" = "Mixed normal",
  "skewed_t" = "Skewed t",
  "contaminated" = "Contaminated",
  "golm" = "GOLM",
  "pareto" = "Pareto",
  "gpd" = "GPD"
)

method_metric_levels <- c(
  "overlap_mis", "overlap_cooks", "overlap_lev", "overlap_dfbetas"
)
method_levels <- c("MIS", "Cook's D", "Leverage", "DFBETAS")

coverage_metric_levels <- c(
  "cover_evd", "cover_cooks", "cover_lev", "cover_dfbetas"
)
coverage_method_levels <- c("MIS-EVT", "Cook's D", "Leverage", "DFBETAS")

contam_values <- NOMINAL_CONTAM_LEVELS
contam_labels <- format_contamination(contam_values)

sim <- sim %>%
  mutate(
    outlier_method = factor(outlier_method, levels = outlier_levels),
    outlier_label = factor(
      unname(outlier_labels[as.character(outlier_method)]),
      levels = unname(outlier_labels[outlier_levels])
    ),
    error_type = factor(error_type, levels = error_levels),
    error_label = factor(
      unname(error_labels[as.character(error_type)]),
      levels = unname(error_labels[error_levels])
    ),
    contam_target = map_to_nominal_contamination(contam_prop),
    contam_label = factor(
      format_contamination(contam_target),
      levels = contam_labels
    )
  )

contamination_mapping_audit <- sim %>%
  distinct(n_obs, set_size, contam_prop, contam_target, contam_label) %>%
  arrange(contam_target, n_obs)

utils::write.csv(
  contamination_mapping_audit,
  file.path(data_dir, "02_contamination_mapping_audit.csv"),
  row.names = FALSE,
  na = ""
)

diag_02b <- diag_02b %>%
  mutate(
    error_type = factor(error_type, levels = error_levels),
    error_label = factor(
      unname(error_labels[as.character(error_type)]),
      levels = unname(error_labels[error_levels])
    ),
    dfb_ratio = safe_ratio(abs(dfb_detected), abs(dfb_injected)),
    leverage_ratio = safe_ratio(mean_lev_detected, mean_lev_injected),
    residual_ratio = safe_ratio(mean_res_detected, mean_res_injected),
    log2_dfb_ratio = safe_log2_ratio(dfb_detected, dfb_injected, use_absolute = TRUE),
    log2_leverage_ratio = safe_log2_ratio(mean_lev_detected, mean_lev_injected),
    log2_residual_ratio = safe_log2_ratio(mean_res_detected, mean_res_injected)
  )

# Iteration-level means within each recorded design cell. This keeps each design
# cell as the unit for distributional summaries rather than treating every row
# as an unrelated observation.
scenario_keys <- intersect(
  c(
    "n_obs", "contam_prop", "set_size", "x_type", "error_type",
    "dist_param", "mix_prop", "outlier_method"
  ),
  names(sim)
)

det_cell <- sim %>%
  filter(as.character(outlier_method) != "none") %>%
  group_by(across(all_of(scenario_keys))) %>%
  summarise(
    overlap_mis = safe_mean(overlap_mis),
    overlap_cooks = safe_mean(overlap_cooks),
    overlap_lev = safe_mean(overlap_lev),
    overlap_dfbetas = safe_mean(overlap_dfbetas),
    detect_mis = safe_mean(detection_success),
    detect_cooks = safe_mean(detect_cooks),
    detect_lev = safe_mean(detect_lev),
    detect_dfbetas = safe_mean(detect_dfbetas),
    power_mis = safe_mean(cover_evd),
    power_cooks = safe_mean(cover_cooks),
    power_lev = safe_mean(cover_lev),
    power_dfbetas = safe_mean(cover_dfbetas),
    convergence = safe_mean(converged),
    mean_compute_time = safe_mean(compute_time),
    n_iterations = n(),
    .groups = "drop"
  ) %>%
  mutate(
    best_classical_overlap = row_max_na(overlap_cooks, overlap_lev, overlap_dfbetas),
    best_classical_power = row_max_na(power_cooks, power_lev, power_dfbetas),
    mis_advantage = overlap_mis - best_classical_overlap,
    outlier_label = factor(
      unname(outlier_labels[as.character(outlier_method)]),
      levels = unname(outlier_labels[outlier_levels[-1L]])
    ),
    error_label = factor(
      unname(error_labels[as.character(error_type)]),
      levels = unname(error_labels[error_levels])
    ),
    contam_target = map_to_nominal_contamination(contam_prop),
    contam_label = factor(
      format_contamination(contam_target),
      levels = contam_labels
    )
  )

utils::write.csv(
  det_cell,
  file.path(data_dir, "02_design_cell_detection_summary.csv"),
  row.names = FALSE,
  na = ""
)

# ==============================================================================
# 4. Main-paper figures
# ==============================================================================

# ----------------------------------------------------------------------------
# Figure 1: MIS advantage over the best classical diagnostic
# ----------------------------------------------------------------------------

heatmap_main <- det_cell %>%
  # Aggregate by the nominal contamination level used for display. The saved
  # realized contam_prop can contain more than one numeric value that maps to
  # the same nominal label, especially when integer k rounding or old cached
  # simulation chunks are present.
  group_by(n_obs, contam_target, contam_label, outlier_label) %>%
  summarise(
    mis_overlap = safe_mean(overlap_mis),
    cooks_overlap = safe_mean(overlap_cooks),
    leverage_overlap = safe_mean(overlap_lev),
    dfbetas_overlap = safe_mean(overlap_dfbetas),
    .groups = "drop"
  ) %>%
  mutate(
    best_classical_overlap = row_max_na(
      cooks_overlap, leverage_overlap, dfbetas_overlap
    ),
    mis_advantage = mis_overlap - best_classical_overlap,
    n_label = factor(n_obs, levels = sort(unique(n_obs)))
  )

max_heat <- max(abs(heatmap_main$mis_advantage), na.rm = TRUE)
if (!is.finite(max_heat) || max_heat == 0) max_heat <- 0.01
heatmap_main <- heatmap_main %>%
  mutate(
    cell_label = sprintf("%+.1f", 100 * mis_advantage),
    text_color = ifelse(abs(mis_advantage) >= 0.58 * max_heat, "white", "black")
  )

utils::write.csv(
  heatmap_main,
  file.path(data_dir, "02_fig1_mis_advantage_heatmap_data.csv"),
  row.names = FALSE,
  na = ""
)

p_fig1 <- ggplot(
  heatmap_main,
  aes(x = n_label, y = contam_label, fill = mis_advantage)
) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = cell_label, colour = text_color), size = 3.25) +
  facet_wrap(~outlier_label, nrow = 1, drop = TRUE) +
  scale_fill_gradient2(
    low = COL_ORANGE,
    mid = COL_NEUTRAL,
    high = COL_BLUE,
    midpoint = 0,
    limits = c(-max_heat, max_heat),
    oob = scales::squish,
    labels = scales::label_number(accuracy = 1, scale = 100, suffix = " pp"),
    name = "MIS advantage over best\nclassical diagnostic"
  ) +
  scale_colour_identity() +
  labs(
    x = "Sample size",
    y = "Contamination proportion"
  ) +
  guides(fill = guide_colourbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = grid::unit(8.0, "cm")
  )) +
  theme_heatmap(base_size = 10.5)

save_plot(
  p_fig1,
  file.path(fig_main_dir, "02_fig1_mis_advantage_heatmap.pdf"),
  width = 9.2,
  height = 4.6
)

# ----------------------------------------------------------------------------
# Figure 2: Distributional robustness of MIS versus best classical diagnostic
# Bars show the 10th-90th percentile across recorded design cells.
# ----------------------------------------------------------------------------

distribution_cell <- det_cell %>%
  select(
    all_of(scenario_keys), outlier_label, error_label,
    overlap_mis, best_classical_overlap
  ) %>%
  pivot_longer(
    cols = c(overlap_mis, best_classical_overlap),
    names_to = "method_key",
    values_to = "overlap"
  ) %>%
  mutate(
    method = factor(
      method_key,
      levels = c("overlap_mis", "best_classical_overlap"),
      labels = c("MIS", "Best classical")
    )
  )

distribution_summary <- distribution_cell %>%
  group_by(error_label, outlier_label, method) %>%
  summarise(
    mean_overlap = safe_mean(overlap),
    q10 = safe_quantile(overlap, 0.10),
    q90 = safe_quantile(overlap, 0.90),
    n_cells = sum(is.finite(overlap)),
    .groups = "drop"
  )

utils::write.csv(
  distribution_summary,
  file.path(data_dir, "02_fig2_distributional_robustness_data.csv"),
  row.names = FALSE,
  na = ""
)

position_methods <- position_dodge(width = 0.62)
p_fig2 <- ggplot(
  distribution_summary,
  aes(x = error_label, y = mean_overlap, colour = method, shape = method)
) +
  geom_errorbar(
    aes(ymin = q10, ymax = q90),
    width = 0.18,
    linewidth = 0.55,
    position = position_methods
  ) +
  geom_point(size = 2.35, stroke = 0.75, position = position_methods) +
  facet_wrap(~outlier_label, nrow = 1, drop = TRUE) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_colour_manual(
    values = c("MIS" = COL_BLUE, "Best classical" = COL_ORANGE),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("MIS" = 16, "Best classical" = 17),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = "MIS top-k recovery with injected set"
  ) +
  theme_paper(base_size = 10.5) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  )

save_plot(
  p_fig2,
  file.path(fig_main_dir, "02_fig2_distributional_robustness.pdf"),
  width = 10.5,
  height = 6.2
)

# ----------------------------------------------------------------------------
# Figure 3: Good-leverage mechanism diagnostic
# Ratio panels use log2 ratios, so zero means equality between the MIS-selected
# and injected sets. The overlap panel uses the original 0-1 scale.
# ----------------------------------------------------------------------------

diag_long <- diag_02b %>%
  transmute(
    error_label,
    `MIS top-k recovery` = overlap_mis,
    `log2 absolute DFBETA ratio` = log2_dfb_ratio,
    `log2 leverage ratio` = log2_leverage_ratio,
    `log2 residual ratio` = log2_residual_ratio
  ) %>%
  pivot_longer(
    cols = -error_label,
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "MIS top-k recovery",
        "log2 absolute DFBETA ratio",
        "log2 leverage ratio",
        "log2 residual ratio"
      ),
      labels = c(
        "MIS top-k recovery",
        "log2 |DFBETA detected / injected|",
        "log2 leverage detected / injected",
        "log2 residual detected / injected"
      )
    )
  )

diag_summary <- diag_long %>%
  group_by(error_label, metric) %>%
  summarise(
    median = safe_median(value),
    q10 = safe_quantile(value, 0.10),
    q90 = safe_quantile(value, 0.90),
    n = sum(is.finite(value)),
    .groups = "drop"
  )

reference_lines <- data.frame(
  metric = factor(
    c(
      "MIS top-k recovery",
      "log2 |DFBETA detected / injected|",
      "log2 leverage detected / injected",
      "log2 residual detected / injected"
    ),
    levels = levels(diag_summary$metric)
  ),
  reference = c(0.90, 0, 0, 0)
)

utils::write.csv(
  diag_summary,
  file.path(data_dir, "02_fig3_good_leverage_mechanism_data.csv"),
  row.names = FALSE,
  na = ""
)

p_fig3 <- ggplot(
  diag_summary,
  aes(x = error_label, y = median)
) +
  geom_hline(
    data = reference_lines,
    aes(yintercept = reference),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.55,
    colour = COL_GREY
  ) +
  geom_errorbar(
    aes(ymin = q10, ymax = q90),
    width = 0.18,
    linewidth = 0.55,
    colour = COL_BLUE
  ) +
  geom_point(size = 2.15, colour = COL_BLUE) +
  facet_wrap(~metric, ncol = 2, scales = "free_y") +
  coord_flip() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_paper(base_size = 10.5) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank()
  )

save_plot(
  p_fig3,
  file.path(fig_main_dir, "02_fig3_good_leverage_mechanism.pdf"),
  width = 9.0,
  height = 7.2
)

# ==============================================================================
# 5. Supplementary figures
# ==============================================================================

# ----------------------------------------------------------------------------
# Figure A1: Absolute overlap heatmaps for all four diagnostics
# ----------------------------------------------------------------------------

overlap_heatmap <- det_cell %>%
  select(
    n_obs, contam_target, contam_label, outlier_label,
    all_of(method_metric_levels)
  ) %>%
  pivot_longer(
    cols = all_of(method_metric_levels),
    names_to = "metric",
    values_to = "overlap"
  ) %>%
  mutate(
    method = factor(metric, levels = method_metric_levels, labels = method_levels)
  ) %>%
  group_by(n_obs, contam_target, contam_label, outlier_label, method) %>%
  summarise(overlap = safe_mean(overlap), .groups = "drop") %>%
  mutate(
    n_label = factor(n_obs, levels = sort(unique(n_obs))),
    cell_label = scales::percent(overlap, accuracy = 1),
    text_color = ifelse(overlap >= 0.62, "white", "black")
  )

utils::write.csv(
  overlap_heatmap,
  file.path(data_dir, "02_figA1_all_method_heatmaps_data.csv"),
  row.names = FALSE,
  na = ""
)

p_figA1 <- ggplot(
  overlap_heatmap,
  aes(x = n_label, y = contam_label, fill = overlap)
) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = cell_label, colour = text_color), size = 2.55) +
  facet_grid(rows = vars(outlier_label), cols = vars(method), drop = TRUE) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = COL_BLUE_DARK,
    limits = c(0, 1),
    oob = scales::squish,
    labels = scales::label_percent(accuracy = 1),
    name = "MIS top-k recovery"
  ) +
  scale_colour_identity() +
  labs(
    x = "Sample size",
    y = "Contamination proportion"
  ) +
  guides(fill = guide_colourbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = grid::unit(8.0, "cm")
  )) +
  theme_heatmap(base_size = 9.5)

save_plot(
  p_figA1,
  file.path(fig_supp_dir, "02_figA1_all_method_overlap_heatmaps.pdf"),
  width = 12.2,
  height = 8.8
)

# ----------------------------------------------------------------------------
# Figure A2: EVT convergence heatmap
# ----------------------------------------------------------------------------

convergence_heatmap <- sim %>%
  group_by(n_obs, contam_target, contam_label, outlier_label) %>%
  summarise(
    convergence = safe_mean(converged),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    n_label = factor(n_obs, levels = sort(unique(n_obs))),
    cell_label = scales::percent(convergence, accuracy = 0.1),
    text_color = ifelse(convergence >= 0.62, "white", "black")
  )

utils::write.csv(
  convergence_heatmap,
  file.path(data_dir, "02_figA2_evt_convergence_heatmap_data.csv"),
  row.names = FALSE,
  na = ""
)

p_figA2 <- ggplot(
  convergence_heatmap,
  aes(x = n_label, y = contam_label, fill = convergence)
) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = cell_label, colour = text_color), size = 3.0) +
  facet_wrap(~outlier_label, ncol = 2, drop = TRUE) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = COL_BLUE_DARK,
    limits = c(0, 1),
    oob = scales::squish,
    labels = scales::label_percent(accuracy = 1),
    name = "EVT convergence"
  ) +
  scale_colour_identity() +
  labs(
    x = "Sample size",
    y = "Contamination proportion"
  ) +
  guides(fill = guide_colourbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = grid::unit(8.0, "cm")
  )) +
  theme_heatmap(base_size = 10.0)

save_plot(
  p_figA2,
  file.path(fig_supp_dir, "02_figA2_evt_convergence_heatmap.pdf"),
  width = 8.6,
  height = 7.0
)

# ----------------------------------------------------------------------------
# Figure A3: Leverage-residual mechanism at the iteration level
# ----------------------------------------------------------------------------

diag_scatter <- diag_02b %>%
  filter(
    is.finite(log2_leverage_ratio),
    is.finite(log2_residual_ratio),
    is.finite(overlap_mis)
  )

utils::write.csv(
  diag_scatter,
  file.path(data_dir, "02_figA3_leverage_residual_scatter_data.csv"),
  row.names = FALSE,
  na = ""
)

p_figA3 <- ggplot(
  diag_scatter,
  aes(
    x = log2_leverage_ratio,
    y = log2_residual_ratio,
    colour = overlap_mis
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45, colour = COL_GREY) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, colour = COL_GREY) +
  geom_point(alpha = 0.58, size = 1.3) +
  facet_wrap(~error_label, ncol = 4, drop = TRUE) +
  scale_colour_gradient(
    low = "#F7FBFF",
    high = COL_BLUE_DARK,
    limits = c(0, 1),
    oob = scales::squish,
    labels = scales::label_percent(accuracy = 1),
    name = "Injected observations recovered\namong MIS top-k selections",
    guide = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(5, "cm"),
      barheight = grid::unit(0.35, "cm")
    )
  ) +
  labs(
    x = "log2 leverage: MIS-selected / injected",
    y = "log2 residual magnitude: MIS-selected / injected"
  ) +
  theme_paper(base_size = 9.5) +
  theme(panel.grid = element_blank())

save_plot(
  p_figA3,
  file.path(fig_supp_dir, "02_figA3_leverage_residual_scatter.pdf"),
  width = 11.0,
  height = 5.8
)

# ==============================================================================
# 6. Main-paper tables
# ==============================================================================

# ----------------------------------------------------------------------------
# Table 1: Detection and rejection summary under contamination
# ----------------------------------------------------------------------------

tab1_raw <- sim %>%
  filter(as.character(outlier_method) != "none") %>%
  group_by(outlier_label) %>%
  summarise(
    mis_overlap = safe_mean(overlap_mis),
    cooks_overlap = safe_mean(overlap_cooks),
    leverage_overlap = safe_mean(overlap_lev),
    dfbetas_overlap = safe_mean(overlap_dfbetas),
    mis_power = safe_mean(cover_evd),
    cooks_power = safe_mean(cover_cooks),
    leverage_power = safe_mean(cover_lev),
    dfbetas_power = safe_mean(cover_dfbetas),
    convergence = safe_mean(converged),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    best_classical_overlap = row_max_na(
      cooks_overlap, leverage_overlap, dfbetas_overlap
    ),
    best_classical_power = row_max_na(
      cooks_power, leverage_power, dfbetas_power
    ),
    mis_advantage = mis_overlap - best_classical_overlap
  )

utils::write.csv(
  tab1_raw,
  file.path(data_dir, "02_tab1_detection_power_summary_numeric.csv"),
  row.names = FALSE,
  na = ""
)

tab1_display <- tab1_raw %>%
  transmute(
    `Outlier mechanism` = as.character(outlier_label),
    `MIS overlap` = fmt_pct(mis_overlap),
    `Best classical overlap` = fmt_pct(best_classical_overlap),
    `MIS advantage` = fmt_pp(mis_advantage),
    `MIS rejection` = fmt_pct(mis_power),
    `Best classical rejection` = fmt_pct(best_classical_power),
    `EVT convergence` = fmt_pct(convergence)
  )

write_tex_table(
  tab1_display,
  tex_path = file.path(tab_main_dir, "02_tab1_detection_power_summary.tex"),
  caption = paste0(
    "Detection overlap and rejection rates across contaminated scenarios. ",
    "The best classical value is the largest aggregate value among Cook's D, ",
    "leverage, and DFBETAS."
  ),
  label = "tab:02-detection-power-summary",
  resize_width = TABLE_WIDTH_WIDE,
  align = "lrrrrrr"
)

# ----------------------------------------------------------------------------
# Table 2: Empirical size under no contamination
# ----------------------------------------------------------------------------

size_long <- sim %>%
  filter(as.character(outlier_method) == "none") %>%
  select(all_of(coverage_metric_levels)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "metric",
    values_to = "reject"
  ) %>%
  mutate(
    method = factor(
      metric,
      levels = coverage_metric_levels,
      labels = coverage_method_levels
    )
  )

size_summary <- size_long %>%
  group_by(method) %>%
  summarise(
    empirical_size = safe_mean(reject),
    deviation_from_05 = empirical_size - 0.05,
    usable_draws = sum(!is.na(reject)),
    .groups = "drop"
  )

utils::write.csv(
  size_summary,
  file.path(data_dir, "02_tab2_empirical_size_summary_numeric.csv"),
  row.names = FALSE,
  na = ""
)

tab2_display <- size_summary %>%
  transmute(
    Method = as.character(method),
    `Empirical rejection` = fmt_pct(empirical_size),
    `Difference from 5%` = fmt_pp(deviation_from_05),
    `Usable draws` = format(usable_draws, big.mark = ",", scientific = FALSE)
  )

write_tex_table(
  tab2_display,
  tex_path = file.path(tab_main_dir, "02_tab2_empirical_size_summary.tex"),
  caption = "Empirical rejection rates under no contamination.",
  label = "tab:02-empirical-size-summary",
  resize_width = TABLE_WIDTH_COMPACT,
  align = "lrrr"
)

# ----------------------------------------------------------------------------
# Table 3: Good-leverage mechanism by error distribution
# ----------------------------------------------------------------------------

diag_table_raw <- diag_02b %>%
  group_by(error_label) %>%
  summarise(
    overlap_median = safe_median(overlap_mis),
    overlap_q10 = safe_quantile(overlap_mis, 0.10),
    overlap_q90 = safe_quantile(overlap_mis, 0.90),
    dfb_ratio_median = safe_median(dfb_ratio),
    dfb_ratio_q10 = safe_quantile(dfb_ratio, 0.10),
    dfb_ratio_q90 = safe_quantile(dfb_ratio, 0.90),
    leverage_ratio_median = safe_median(leverage_ratio),
    leverage_ratio_q10 = safe_quantile(leverage_ratio, 0.10),
    leverage_ratio_q90 = safe_quantile(leverage_ratio, 0.90),
    residual_ratio_median = safe_median(residual_ratio),
    residual_ratio_q10 = safe_quantile(residual_ratio, 0.10),
    residual_ratio_q90 = safe_quantile(residual_ratio, 0.90),
    n = n(),
    .groups = "drop"
  )

utils::write.csv(
  diag_table_raw,
  file.path(data_dir, "02_tab3_good_leverage_mechanism_numeric.csv"),
  row.names = FALSE,
  na = ""
)

tab3_display <- diag_table_raw %>%
  transmute(
    `Error distribution` = as.character(error_label),
    `MIS overlap` = fmt_pct_interval(
      overlap_median, overlap_q10, overlap_q90, digits = 1
    ),
    `Absolute DFBETA ratio` = fmt_interval(
      dfb_ratio_median, dfb_ratio_q10, dfb_ratio_q90, digits = 2
    ),
    `Leverage ratio` = fmt_interval(
      leverage_ratio_median, leverage_ratio_q10, leverage_ratio_q90, digits = 2
    ),
    `Residual ratio` = fmt_interval(
      residual_ratio_median, residual_ratio_q10, residual_ratio_q90, digits = 2
    )
  )

write_tex_table(
  tab3_display,
  tex_path = file.path(tab_main_dir, "02_tab3_good_leverage_mechanism.tex"),
  caption = paste0(
    "Good-leverage mechanism diagnostic by error distribution. Entries are ",
    "medians with 10th--90th percentile intervals; ratios compare the ",
    "MIS-selected set with the injected set."
  ),
  label = "tab:02-good-leverage-mechanism",
  resize_width = TABLE_WIDTH_WIDE,
  align = "lrrrr"
)

# ==============================================================================
# 7. Supplementary tables
# ==============================================================================

# ----------------------------------------------------------------------------
# Table A1: MIS advantage grid by sample size and contamination proportion
# ----------------------------------------------------------------------------

heatmap_key_audit <- heatmap_main %>%
  count(n_obs, contam_label, outlier_label, name = "n_rows") %>%
  filter(n_rows > 1L)

utils::write.csv(
  heatmap_key_audit,
  file.path(diag_dir, "02_heatmap_duplicate_key_audit.csv"),
  row.names = FALSE,
  na = ""
)

tabA1_raw <- heatmap_main %>%
  group_by(n_obs, contam_label, outlier_label) %>%
  summarise(
    mis_advantage = safe_mean(mis_advantage),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = outlier_label,
    values_from = mis_advantage
  ) %>%
  arrange(n_obs, contam_label)

utils::write.csv(
  tabA1_raw,
  file.path(data_dir, "02_tabA1_mis_advantage_grid_numeric.csv"),
  row.names = FALSE,
  na = ""
)

tabA1_value_columns <- setdiff(names(tabA1_raw), c("n_obs", "contam_label"))
tabA1_display <- tabA1_raw %>%
  mutate(across(all_of(tabA1_value_columns), ~fmt_pp(.x))) %>%
  rename(
    `Sample size` = n_obs,
    `Contamination` = contam_label
  )

write_tex_table(
  tabA1_display,
  tex_path = file.path(tab_supp_dir, "02_tabA1_mis_advantage_grid.tex"),
  caption = paste0(
    "MIS overlap advantage, in percentage points, over the best classical ",
    "diagnostic by sample size and contamination proportion."
  ),
  label = "tab:02A-mis-advantage-grid",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = paste0("ll", paste(rep("r", ncol(tabA1_display) - 2L), collapse = ""))
)

# ----------------------------------------------------------------------------
# Table A2: EVT convergence grid
# ----------------------------------------------------------------------------

tabA2_raw <- convergence_heatmap %>%
  group_by(n_obs, contam_label, outlier_label) %>%
  summarise(
    convergence = safe_mean(convergence),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = outlier_label,
    values_from = convergence
  ) %>%
  arrange(n_obs, contam_label)

utils::write.csv(
  tabA2_raw,
  file.path(data_dir, "02_tabA2_evt_convergence_grid_numeric.csv"),
  row.names = FALSE,
  na = ""
)

tabA2_value_columns <- setdiff(names(tabA2_raw), c("n_obs", "contam_label"))
tabA2_display <- tabA2_raw %>%
  mutate(across(all_of(tabA2_value_columns), ~fmt_pct(.x))) %>%
  rename(
    `Sample size` = n_obs,
    `Contamination` = contam_label
  )

write_tex_table(
  tabA2_display,
  tex_path = file.path(tab_supp_dir, "02_tabA2_evt_convergence_grid.tex"),
  caption = "EVT convergence rates by sample size, contamination proportion, and contamination mechanism.",
  label = "tab:02A-evt-convergence-grid",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = paste0("ll", paste(rep("r", ncol(tabA2_display) - 2L), collapse = ""))
)

# ----------------------------------------------------------------------------
# Table A3: Error-distribution-specific detection summary
# ----------------------------------------------------------------------------

tabA3_raw <- sim %>%
  filter(as.character(outlier_method) != "none") %>%
  group_by(error_label, outlier_label) %>%
  summarise(
    mis_overlap = safe_mean(overlap_mis),
    cooks_overlap = safe_mean(overlap_cooks),
    leverage_overlap = safe_mean(overlap_lev),
    dfbetas_overlap = safe_mean(overlap_dfbetas),
    convergence = safe_mean(converged),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    best_classical_overlap = row_max_na(
      cooks_overlap, leverage_overlap, dfbetas_overlap
    ),
    mis_advantage = mis_overlap - best_classical_overlap
  )

utils::write.csv(
  tabA3_raw,
  file.path(data_dir, "02_tabA3_detection_by_error_numeric.csv"),
  row.names = FALSE,
  na = ""
)

tabA3_display <- tabA3_raw %>%
  transmute(
    `Error distribution` = as.character(error_label),
    `Outlier mechanism` = as.character(outlier_label),
    `MIS overlap` = fmt_pct(mis_overlap),
    `Best classical overlap` = fmt_pct(best_classical_overlap),
    `MIS advantage` = fmt_pp(mis_advantage),
    `EVT convergence` = fmt_pct(convergence)
  )

write_tex_table(
  tabA3_display,
  tex_path = file.path(tab_supp_dir, "02_tabA3_detection_by_error.tex"),
  caption = "Detection overlap and EVT convergence by error distribution and contamination mechanism.",
  label = "tab:02A-detection-by-error",
  resize_width = TABLE_WIDTH_WIDE,
  align = "llrrrr"
)

# ----------------------------------------------------------------------------
# Table A4: GEV parameter and rejection summary by error distribution
# ----------------------------------------------------------------------------

tabA4_raw <- sim %>%
  group_by(error_label, outlier_label) %>%
  summarise(
    convergence = safe_mean(converged),
    rejection = safe_mean(cover_evd),
    shape_median = safe_median(shape[converged %in% TRUE]),
    scale_median = safe_median(scale[converged %in% TRUE]),
    location_median = safe_median(loc[converged %in% TRUE]),
    median_compute_time = safe_median(compute_time),
    n = n(),
    .groups = "drop"
  )

utils::write.csv(
  tabA4_raw,
  file.path(data_dir, "02_tabA4_gev_summary_numeric.csv"),
  row.names = FALSE,
  na = ""
)

tabA4_display <- tabA4_raw %>%
  transmute(
    `Error distribution` = as.character(error_label),
    `Scenario` = as.character(outlier_label),
    Convergence = fmt_pct(convergence),
    Rejection = fmt_pct(rejection),
    `Median shape` = fmt_num(shape_median, 3),
    `Median scale` = fmt_num(scale_median, 3),
    `Median location` = fmt_num(location_median, 3),
    `Median time (s)` = fmt_num(median_compute_time, 3)
  )

write_tex_table(
  tabA4_display,
  tex_path = file.path(tab_supp_dir, "02_tabA4_gev_summary.tex"),
  caption = "GEV convergence, rejection, parameter, and computation-time summaries by error distribution and scenario.",
  label = "tab:02A-gev-summary",
  resize_width = TABLE_WIDTH_WIDE,
  align = "llrrrrrr"
)

# ----------------------------------------------------------------------------
# Table A5: Adaptive block-count summary
# ----------------------------------------------------------------------------

tabA5_raw <- sim %>%
  group_by(n_obs, contam_prop, contam_label, set_size) %>%
  summarise(
    block_min = min(block_count, na.rm = TRUE),
    block_median = safe_median(block_count),
    block_max = max(block_count, na.rm = TRUE),
    low_or_infeasible_rate = safe_mean(block_count < 10),
    .groups = "drop"
  ) %>%
  arrange(n_obs, contam_prop)

utils::write.csv(
  tabA5_raw,
  file.path(data_dir, "02_tabA5_block_count_summary_numeric.csv"),
  row.names = FALSE,
  na = ""
)

tabA5_display <- tabA5_raw %>%
  transmute(
    `Sample size` = n_obs,
    `Contamination` = as.character(contam_label),
    `Set size` = set_size,
    `Minimum blocks` = block_min,
    `Median blocks` = fmt_num(block_median, 0),
    `Maximum blocks` = block_max,
    `Blocks below 10` = fmt_pct(low_or_infeasible_rate)
  )

write_tex_table(
  tabA5_display,
  tex_path = file.path(tab_supp_dir, "02_tabA5_block_count_summary.tex"),
  caption = "Adaptive block-count summary by sample size and contamination proportion.",
  label = "tab:02A-block-count-summary",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrrrr"
)

# Full design-cell data are retained in CSV because a complete LaTeX rendering
# would be too tall for a single formal-paper table environment.
utils::write.csv(
  det_cell,
  file.path(tab_supp_dir, "02_tabA6_full_design_cell_results.csv"),
  row.names = FALSE,
  na = ""
)

# ==============================================================================
# 8. Diagnostics and manifest
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(diag_dir, "02_output_session_info.txt")
)

manifest_path <- file.path(diag_dir, "02_output_manifest.csv")
all_outputs <- list.files(
  output_root,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)
all_outputs <- setdiff(normalizePath(all_outputs, winslash = "/", mustWork = FALSE),
                       normalizePath(manifest_path, winslash = "/", mustWork = FALSE))
manifest <- data.frame(
  file = substring(all_outputs, nchar(output_root) + 2L),
  bytes = file.info(all_outputs)$size,
  modified = as.character(file.info(all_outputs)$mtime),
  stringsAsFactors = FALSE
)
utils::write.csv(
  manifest,
  manifest_path,
  row.names = FALSE
)

cat("\nScript 82 completed successfully.\n")
cat("Output directory:\n  ", output_root, "\n", sep = "")
cat(sprintf("Generated %d output files.\n", nrow(manifest) + 1L))
