# ==============================================================================
# File: scripts/83_output_algorithmic.R
# Purpose:
#   Generate publication-ready figures and tables for:
#     1. Script 03v2: algorithmic accuracy, influence capture, and runtime
#     2. Script 03v2: nestedness sub-analysis
#     3. Script 03b: null calibration of constrained and robust GEV fitting
#
# Inputs:
#   output/03v2_scaling_results_master.rds
#   output/03v2_nestedness_traces.rds
#   output/03b_null_calibration.rds
#
# Output structure:
#   output/03_algorithmic/
#     figures/main/          Main-paper PDF figures
#     figures/supplement/    Supplementary PDF figures
#     tables/main/           Complete LaTeX table environments + CSV companions
#     tables/supplement/     Supplementary LaTeX tables + CSV companions
#     data/                  Data used to construct figures and tables
#     diagnostics/           Session information, audits, and output manifest
#
# LaTeX requirements for generated tables:
#   \usepackage{float}
#   \usepackage{graphicx}
#   \usepackage{booktabs}
#
# Design rules:
#   - No title or subtitle inside plots; captions are handled in LaTeX.
#   - Labels, legends, facets, and margins are sized for formal papers.
#   - Color is never the only indicator: line type and point shape are also used.
#   - Diverging heatmaps use orange-white-blue, never red-green.
#   - Main tables are complete \begin{table}[H] environments.
#   - Compact tables use 0.7\columnwidth; wider tables use a larger width only
#     when needed to preserve readable labels and indicators.
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
    ". Install them before running scripts/83_output_algorithmic.R."
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

input_scaling <- file.path(
  project_root, "output", "03v2_scaling_results_master.rds"
)
input_nested <- file.path(
  project_root, "output", "03v2_nestedness_traces.rds"
)
input_null <- file.path(
  project_root, "output", "03b_null_calibration.rds"
)

output_root <- file.path(project_root, "output", "03_algorithmic")
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

if (!file.exists(input_scaling)) {
  stop("Missing input file: ", input_scaling)
}
if (!file.exists(input_null)) {
  stop("Missing input file: ", input_null)
}
if (!file.exists(input_nested)) {
  warning(
    "Nestedness input is missing: ", input_nested,
    ". Nestedness figures and tables will be skipped."
  )
}

# PDF is the publication format. Set TRUE only when raster previews are useful.
SAVE_PNG_PREVIEWS <- FALSE
PNG_DPI <- 320

# Color-blind-conscious palette. Heatmaps use orange-white-blue.
COL_ORANGE <- "#E69F00"
COL_ORANGE_DARK <- "#D55E00"
COL_BLUE <- "#0072B2"
COL_BLUE_LIGHT <- "#56B4E9"
COL_BLUE_DARK <- "#08519C"
COL_BLACK <- "#111111"
COL_GREY <- "#7A7A7A"
COL_GREY_LIGHT <- "#B8B8B8"
COL_NEUTRAL <- "#F7F7F7"

TABLE_WIDTH_COMPACT <- "0.7\\columnwidth"
TABLE_WIDTH_MEDIUM  <- "0.85\\columnwidth"
TABLE_WIDTH_WIDE    <- "\\columnwidth"

# Empirical size is evaluated at these nominal levels.
ALPHA_GRID <- c(0.01, 0.025, 0.05, 0.10)

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
  unname(stats::quantile(
    x, probs = probability, na.rm = TRUE, names = FALSE, type = 7
  ))
}

safe_ks_pvalue <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 20L) return(NA_real_)
  suppressWarnings(stats::ks.test(x, "punif")$p.value)
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
    placement = "H",
    font_command = NULL
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
    if (!is.null(font_command)) font_command else character(0),
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
  pdf_path <- file.path(
    dirname(filename),
    paste0(tools::file_path_sans_ext(basename(filename)), ".pdf")
  )
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
    png_path <- file.path(
      dirname(filename),
      paste0(tools::file_path_sans_ext(basename(filename)), ".png")
    )
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
      legend.key.width = grid::unit(1.25, "cm"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      panel.spacing = grid::unit(1.25, "lines"),
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

scaling <- readRDS(input_scaling)
null_cal <- readRDS(input_null)
nested <- if (file.exists(input_nested)) readRDS(input_nested) else data.frame()

required_scaling_columns <- c(
  "iter", "N", "k", "B_type", "B_actual", "B_raw", "B_capped",
  "architecture", "rho",
  "det_greedy_noref", "det_greedy_ref", "det_exact_noref", "det_exact_ref",
  "ir_greedy_noref", "ir_greedy_ref", "ir_exact_noref", "ir_exact_ref",
  "cpu_det_greedy_noref", "cpu_det_greedy_ref",
  "cpu_det_exact_noref", "cpu_det_exact_ref",
  "p_greedy_con", "conv_greedy_con", "shape_greedy_con",
  "p_exact_con", "conv_exact_con", "shape_exact_con",
  "p_exact_rob", "conv_exact_rob", "shape_exact_rob",
  "cpu_evt_greedy_con", "cpu_evt_exact_con", "cpu_evt_exact_rob",
  "error_msg"
)
missing_scaling_columns <- setdiff(required_scaling_columns, names(scaling))
if (length(missing_scaling_columns) > 0L) {
  stop(
    "03v2_scaling_results_master.rds is missing required column(s): ",
    paste(missing_scaling_columns, collapse = ", ")
  )
}

required_null_columns <- c(
  "iter", "N", "k", "B", "architecture",
  "p_constrained", "conv_constrained", "shape_constrained",
  "p_robust", "conv_robust", "shape_robust", "error_msg"
)
missing_null_columns <- setdiff(required_null_columns, names(null_cal))
if (length(missing_null_columns) > 0L) {
  stop(
    "03b_null_calibration.rds is missing required column(s): ",
    paste(missing_null_columns, collapse = ", ")
  )
}

if (nrow(nested) > 0L) {
  required_nested_columns <- c(
    "k", "jaccard_with_prev", "nested", "influence_magnitude",
    "cpu_seconds", "architecture", "method"
  )
  missing_nested_columns <- setdiff(required_nested_columns, names(nested))
  if (length(missing_nested_columns) > 0L) {
    warning(
      "03v2_nestedness_traces.rds is missing required column(s): ",
      paste(missing_nested_columns, collapse = ", "),
      ". Nestedness outputs will be skipped."
    )
    nested <- data.frame()
  }
}

cat(sprintf("Loaded Script 03v2 scaling data: %s rows\n",
            format(nrow(scaling), big.mark = ",")))
cat(sprintf("Loaded Script 03b null data:     %s rows\n",
            format(nrow(null_cal), big.mark = ",")))
if (nrow(nested) > 0L) {
  cat(sprintf("Loaded nestedness traces:        %s rows\n",
              format(nrow(nested), big.mark = ",")))
}

# ==============================================================================
# 3. Labels and tidy preparation
# ==============================================================================

architecture_order <- c(
  "simple", "complex",
  "interaction", "triple_interaction", "sparse_binary_interaction",
  "polynomial_interaction", "high_k_interaction",
  "nonlinear_nuisance", "collinear_interaction",
  "plm_confounded", "plm_nonlinear"
)

architecture_labels <- c(
  "simple" = "Simple linear",
  "complex" = "Complex linear",
  "interaction" = "Interaction",
  "triple_interaction" = "Triple interaction",
  "sparse_binary_interaction" = "Sparse binary interaction",
  "polynomial_interaction" = "Polynomial interaction",
  "high_k_interaction" = "High-k interaction",
  "nonlinear_nuisance" = "Nonlinear nuisance",
  "collinear_interaction" = "Collinear interaction",
  "plm_confounded" = "PLM: confounded",
  "plm_nonlinear" = "PLM: nonlinear"
)

architecture_family <- c(
  "simple" = "Linear designs",
  "complex" = "Linear designs",
  "interaction" = "Interaction designs",
  "triple_interaction" = "Interaction designs",
  "sparse_binary_interaction" = "Interaction designs",
  "polynomial_interaction" = "Interaction designs",
  "high_k_interaction" = "Interaction designs",
  "nonlinear_nuisance" = "Difficult designs",
  "collinear_interaction" = "Difficult designs",
  "plm_confounded" = "Partially linear designs",
  "plm_nonlinear" = "Partially linear designs"
)

family_order <- c(
  "Linear designs", "Interaction designs",
  "Difficult designs", "Partially linear designs"
)

det_method_order <- c(
  "greedy_noref", "greedy_ref", "exact_noref", "exact_ref"
)
det_method_labels <- c(
  "greedy_noref" = "Greedy",
  "greedy_ref" = "Greedy + refinement",
  "exact_noref" = "Exact",
  "exact_ref" = "Exact + refinement"
)
det_method_colors <- c(
  "Greedy" = COL_GREY,
  "Greedy + refinement" = COL_ORANGE,
  "Exact" = COL_BLACK,
  "Exact + refinement" = COL_BLUE
)
det_method_linetypes <- c(
  "Greedy" = "dotted",
  "Greedy + refinement" = "dashed",
  "Exact" = "dotdash",
  "Exact + refinement" = "solid"
)
det_method_shapes <- c(
  "Greedy" = 1,
  "Greedy + refinement" = 17,
  "Exact" = 4,
  "Exact + refinement" = 16
)

evt_method_order <- c("greedy_con", "exact_con", "exact_rob")
evt_method_labels <- c(
  "greedy_con" = "Greedy + constrained GEV",
  "exact_con" = "Exact + constrained GEV",
  "exact_rob" = "Exact + robust GEV"
)
evt_method_colors <- c(
  "Greedy + constrained GEV" = COL_GREY,
  "Exact + constrained GEV" = COL_ORANGE_DARK,
  "Exact + robust GEV" = COL_BLUE
)
evt_method_linetypes <- c(
  "Greedy + constrained GEV" = "dotted",
  "Exact + constrained GEV" = "dashed",
  "Exact + robust GEV" = "solid"
)
evt_method_shapes <- c(
  "Greedy + constrained GEV" = 1,
  "Exact + constrained GEV" = 17,
  "Exact + robust GEV" = 16
)

null_method_order <- c("constrained", "robust")
null_method_labels <- c(
  "constrained" = "Constrained GEV",
  "robust" = "Robust GEV"
)
null_method_colors <- c(
  "Constrained GEV" = COL_ORANGE,
  "Robust GEV" = COL_BLUE
)
null_method_linetypes <- c(
  "Constrained GEV" = "dashed",
  "Robust GEV" = "solid"
)
null_method_shapes <- c(
  "Constrained GEV" = 17,
  "Robust GEV" = 16
)

add_architecture_labels <- function(df) {
  df %>%
    mutate(
      architecture = as.character(architecture),
      architecture_label = ifelse(
        architecture %in% names(architecture_labels),
        unname(architecture_labels[architecture]),
        architecture
      ),
      architecture_label = factor(
        architecture_label,
        levels = unique(c(
          unname(architecture_labels[architecture_order]),
          architecture_label
        ))
      ),
      family_label = ifelse(
        architecture %in% names(architecture_family),
        unname(architecture_family[architecture]),
        "Other designs"
      ),
      family_label = factor(
        family_label,
        levels = c(family_order, "Other designs")
      )
    )
}

scaling <- add_architecture_labels(scaling)
null_cal <- add_architecture_labels(null_cal)
if (nrow(nested) > 0L) nested <- add_architecture_labels(nested)

# Detection data in long form: one row per Monte Carlo iteration and method.
det_spec <- data.frame(
  method_id = det_method_order,
  overlap_col = c(
    "det_greedy_noref", "det_greedy_ref", "det_exact_noref", "det_exact_ref"
  ),
  ir_col = c(
    "ir_greedy_noref", "ir_greedy_ref", "ir_exact_noref", "ir_exact_ref"
  ),
  cpu_col = c(
    "cpu_det_greedy_noref", "cpu_det_greedy_ref",
    "cpu_det_exact_noref", "cpu_det_exact_ref"
  ),
  stringsAsFactors = FALSE
)

det_long <- bind_rows(lapply(seq_len(nrow(det_spec)), function(i) {
  scaling %>%
    transmute(
      iter, N, k, B_type, B_actual, B_raw, B_capped,
      architecture, architecture_label, family_label, rho, error_msg,
      method_id = det_spec$method_id[[i]],
      overlap = .data[[det_spec$overlap_col[[i]]]],
      influence_ratio = .data[[det_spec$ir_col[[i]]]],
      cpu_seconds = .data[[det_spec$cpu_col[[i]]]]
    )
})) %>%
  mutate(
    method_id = factor(method_id, levels = det_method_order),
    method_label = factor(
      unname(det_method_labels[as.character(method_id)]),
      levels = unname(det_method_labels[det_method_order])
    )
  )

# EVT data in long form.
evt_spec <- data.frame(
  method_id = evt_method_order,
  p_col = c("p_greedy_con", "p_exact_con", "p_exact_rob"),
  conv_col = c("conv_greedy_con", "conv_exact_con", "conv_exact_rob"),
  shape_col = c("shape_greedy_con", "shape_exact_con", "shape_exact_rob"),
  cpu_col = c("cpu_evt_greedy_con", "cpu_evt_exact_con", "cpu_evt_exact_rob"),
  stringsAsFactors = FALSE
)

evt_long <- bind_rows(lapply(seq_len(nrow(evt_spec)), function(i) {
  scaling %>%
    transmute(
      iter, N, k, B_type, B_actual, B_raw, B_capped,
      architecture, architecture_label, family_label, rho, error_msg,
      method_id = evt_spec$method_id[[i]],
      p_value = .data[[evt_spec$p_col[[i]]]],
      converged = .data[[evt_spec$conv_col[[i]]]],
      shape = .data[[evt_spec$shape_col[[i]]]],
      cpu_seconds = .data[[evt_spec$cpu_col[[i]]]]
    )
})) %>%
  mutate(
    method_id = factor(method_id, levels = evt_method_order),
    method_label = factor(
      unname(evt_method_labels[as.character(method_id)]),
      levels = unname(evt_method_labels[evt_method_order])
    )
  )

# Null-calibration data in long form.
null_long <- bind_rows(
  null_cal %>%
    transmute(
      iter, N, k, B, architecture, architecture_label, family_label, error_msg,
      method_id = "constrained",
      p_value = p_constrained,
      converged = conv_constrained,
      shape = shape_constrained
    ),
  null_cal %>%
    transmute(
      iter, N, k, B, architecture, architecture_label, family_label, error_msg,
      method_id = "robust",
      p_value = p_robust,
      converged = conv_robust,
      shape = shape_robust
    )
) %>%
  mutate(
    method_id = factor(method_id, levels = null_method_order),
    method_label = factor(
      unname(null_method_labels[as.character(method_id)]),
      levels = unname(null_method_labels[null_method_order])
    )
  )

# Scenario-level summaries. Each design cell, rather than each raw iteration,
# is the unit used when aggregating across architectures and tuning parameters.
det_cell <- det_long %>%
  group_by(
    architecture, architecture_label, family_label,
    N, k, B_type, B_actual, B_raw, B_capped, rho,
    method_id, method_label
  ) %>%
  summarise(
    overlap_mean = safe_mean(overlap),
    overlap_median = safe_median(overlap),
    influence_ratio_median = safe_median(influence_ratio),
    influence_ratio_q10 = safe_quantile(influence_ratio, 0.10),
    influence_ratio_q90 = safe_quantile(influence_ratio, 0.90),
    cpu_median = safe_median(cpu_seconds[cpu_seconds > 0]),
    cpu_q25 = safe_quantile(cpu_seconds[cpu_seconds > 0], 0.25),
    cpu_q75 = safe_quantile(cpu_seconds[cpu_seconds > 0], 0.75),
    n_iterations = n(),
    n_overlap = sum(is.finite(overlap)),
    n_ir = sum(is.finite(influence_ratio)),
    n_cpu = sum(is.finite(cpu_seconds) & cpu_seconds > 0),
    failure_rate = safe_mean(!is.na(error_msg) & nzchar(as.character(error_msg))),
    .groups = "drop"
  )

evt_cell <- evt_long %>%
  group_by(
    architecture, architecture_label, family_label,
    N, k, B_type, B_actual, B_raw, B_capped, rho,
    method_id, method_label
  ) %>%
  summarise(
    convergence_rate = safe_mean(as.numeric(converged)),
    rejection_05 = safe_mean(as.numeric(p_value < 0.05)),
    shape_median = safe_median(shape[converged %in% TRUE]),
    shape_q25 = safe_quantile(shape[converged %in% TRUE], 0.25),
    shape_q75 = safe_quantile(shape[converged %in% TRUE], 0.75),
    cpu_median = safe_median(cpu_seconds[cpu_seconds > 0]),
    cpu_q25 = safe_quantile(cpu_seconds[cpu_seconds > 0], 0.25),
    cpu_q75 = safe_quantile(cpu_seconds[cpu_seconds > 0], 0.75),
    n_iterations = n(),
    n_valid_p = sum(is.finite(p_value)),
    n_converged = sum(converged %in% TRUE, na.rm = TRUE),
    .groups = "drop"
  )

null_cell <- null_long %>%
  group_by(
    architecture, architecture_label, family_label,
    N, k, B, method_id, method_label
  ) %>%
  summarise(
    convergence_rate = safe_mean(as.numeric(converged)),
    shape_median = safe_median(shape[converged %in% TRUE]),
    shape_q25 = safe_quantile(shape[converged %in% TRUE], 0.25),
    shape_q75 = safe_quantile(shape[converged %in% TRUE], 0.75),
    ks_pvalue = safe_ks_pvalue(p_value),
    n_iterations = n(),
    n_valid_p = sum(is.finite(p_value)),
    .groups = "drop"
  )

null_alpha_cell <- bind_rows(lapply(ALPHA_GRID, function(alpha_value) {
  null_long %>%
    group_by(
      architecture, architecture_label, family_label,
      N, k, B, method_id, method_label
    ) %>%
    summarise(
      rejection_rate = safe_mean(as.numeric(p_value < alpha_value)),
      n_valid_p = sum(is.finite(p_value)),
      .groups = "drop"
    ) %>%
    mutate(alpha = alpha_value)
}))

utils::write.csv(
  det_cell,
  file.path(data_dir, "03_design_cell_detection_summary.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  evt_cell,
  file.path(data_dir, "03_design_cell_evt_summary.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  null_cell,
  file.path(data_dir, "03_design_cell_null_summary.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  null_alpha_cell,
  file.path(data_dir, "03_design_cell_null_rejection_rates.csv"),
  row.names = FALSE,
  na = ""
)

# ==============================================================================
# 4. Main-paper figures
# ==============================================================================

# ----------------------------------------------------------------------------
# Figure 1: Detection overlap across set sizes and architecture families
# Bands show the 10th-90th percentile across recorded design cells.
# ----------------------------------------------------------------------------

detection_main <- det_cell %>%
  group_by(family_label, k, method_id, method_label) %>%
  summarise(
    overlap = safe_mean(overlap_mean),
    lower = safe_quantile(overlap_mean, 0.10),
    upper = safe_quantile(overlap_mean, 0.90),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  filter(is.finite(overlap))

utils::write.csv(
  detection_main,
  file.path(data_dir, "03_fig1_detection_overlap_data.csv"),
  row.names = FALSE,
  na = ""
)

p_fig1 <- ggplot(
  detection_main,
  aes(
    x = k, y = overlap,
    colour = method_label, fill = method_label,
    linetype = method_label, shape = method_label,
    group = method_label
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.08, colour = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2.2, stroke = 0.8) +
  facet_wrap(~family_label, ncol = 2, scales = "free_x", drop = TRUE) +
  scale_colour_manual(values = det_method_colors, drop = FALSE) +
  scale_fill_manual(values = det_method_colors, drop = FALSE) +
  scale_linetype_manual(values = det_method_linetypes, drop = FALSE) +
  scale_shape_manual(values = det_method_shapes, drop = FALSE) +
  scale_x_continuous(breaks = sort(unique(detection_main$k))) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  labs(
    x = "Influential-set size, k",
    y = "Mean detection overlap",
    colour = NULL, fill = NULL, linetype = NULL, shape = NULL
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_paper(base_size = 10.5)

save_plot(
  p_fig1,
  file.path(fig_main_dir, "03_fig1_detection_overlap.pdf"),
  width = 8.8,
  height = 6.7
)

# ----------------------------------------------------------------------------
# Figure 2a: Detection runtime scaling
# ----------------------------------------------------------------------------

detection_runtime <- det_cell %>%
  filter(is.finite(cpu_median), cpu_median > 0) %>%
  group_by(N, method_id, method_label) %>%
  summarise(
    median_seconds = safe_median(cpu_median),
    lower = safe_quantile(cpu_median, 0.25),
    upper = safe_quantile(cpu_median, 0.75),
    n_cells = n(),
    .groups = "drop"
  )

utils::write.csv(
  detection_runtime,
  file.path(data_dir, "03_fig2a_detection_runtime_data.csv"),
  row.names = FALSE,
  na = ""
)

p_fig2a <- ggplot(
  detection_runtime,
  aes(
    x = N, y = median_seconds,
    colour = method_label, fill = method_label,
    linetype = method_label, shape = method_label,
    group = method_label
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.08, colour = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.4, stroke = 0.8) +
  scale_colour_manual(values = det_method_colors, drop = FALSE) +
  scale_fill_manual(values = det_method_colors, drop = FALSE) +
  scale_linetype_manual(values = det_method_linetypes, drop = FALSE) +
  scale_shape_manual(values = det_method_shapes, drop = FALSE) +
  scale_x_continuous(
    breaks = sort(unique(detection_runtime$N)),
    labels = scales::label_number(big.mark = ",")
  ) +
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.001),
    expand = expansion(mult = c(0.04, 0.10))
  ) +
  labs(
    x = "Sample size",
    y = "Median detection time (seconds, log scale)",
    colour = NULL, fill = NULL, linetype = NULL, shape = NULL
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_paper(base_size = 10.5)

save_plot(
  p_fig2a,
  file.path(fig_main_dir, "03_fig2a_detection_runtime.pdf"),
  width = 6.8,
  height = 4.7
)

# ----------------------------------------------------------------------------
# Figure 2b: EVT runtime scaling
# ----------------------------------------------------------------------------

evt_runtime <- evt_cell %>%
  filter(is.finite(cpu_median), cpu_median > 0) %>%
  group_by(N, method_id, method_label) %>%
  summarise(
    median_seconds = safe_median(cpu_median),
    lower = safe_quantile(cpu_median, 0.25),
    upper = safe_quantile(cpu_median, 0.75),
    n_cells = n(),
    .groups = "drop"
  )

utils::write.csv(
  evt_runtime,
  file.path(data_dir, "03_fig2b_evt_runtime_data.csv"),
  row.names = FALSE,
  na = ""
)

p_fig2b <- ggplot(
  evt_runtime,
  aes(
    x = N, y = median_seconds,
    colour = method_label, fill = method_label,
    linetype = method_label, shape = method_label,
    group = method_label
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.08, colour = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.4, stroke = 0.8) +
  scale_colour_manual(values = evt_method_colors, drop = FALSE) +
  scale_fill_manual(values = evt_method_colors, drop = FALSE) +
  scale_linetype_manual(values = evt_method_linetypes, drop = FALSE) +
  scale_shape_manual(values = evt_method_shapes, drop = FALSE) +
  scale_x_continuous(
    breaks = sort(unique(evt_runtime$N)),
    labels = scales::label_number(big.mark = ",")
  ) +
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.001),
    expand = expansion(mult = c(0.04, 0.10))
  ) +
  labs(
    x = "Sample size",
    y = "Median EVT time (seconds, log scale)",
    colour = NULL, fill = NULL, linetype = NULL, shape = NULL
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_paper(base_size = 10.5)

save_plot(
  p_fig2b,
  file.path(fig_main_dir, "03_fig2b_evt_runtime.pdf"),
  width = 6.8,
  height = 4.7
)

# ----------------------------------------------------------------------------
# Figure 3: Null calibration curves
# Bands show the 10th-90th percentile across design cells.
# ----------------------------------------------------------------------------

calibration_main <- null_alpha_cell %>%
  group_by(family_label, method_id, method_label, alpha) %>%
  summarise(
    empirical_size = safe_mean(rejection_rate),
    lower = safe_quantile(rejection_rate, 0.10),
    upper = safe_quantile(rejection_rate, 0.90),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  filter(is.finite(empirical_size))

utils::write.csv(
  calibration_main,
  file.path(data_dir, "03_fig3_null_calibration_data.csv"),
  row.names = FALSE,
  na = ""
)

p_fig3 <- ggplot(
  calibration_main,
  aes(
    x = alpha, y = empirical_size,
    colour = method_label, fill = method_label,
    linetype = method_label, shape = method_label,
    group = method_label
  )
) +
  geom_abline(
    intercept = 0, slope = 1,
    colour = COL_GREY, linewidth = 0.55, linetype = "longdash"
  ) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.10, colour = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.4, stroke = 0.8) +
  facet_wrap(~family_label, ncol = 2, drop = TRUE) +
  scale_colour_manual(values = null_method_colors, drop = FALSE) +
  scale_fill_manual(values = null_method_colors, drop = FALSE) +
  scale_linetype_manual(values = null_method_linetypes, drop = FALSE) +
  scale_shape_manual(values = null_method_shapes, drop = FALSE) +
  scale_x_continuous(
    breaks = ALPHA_GRID,
    labels = scales::label_percent(accuracy = 0.1),
    expand = expansion(mult = c(0.04, 0.06))
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = "Nominal significance level",
    y = "Empirical rejection rate",
    colour = NULL, fill = NULL, linetype = NULL, shape = NULL
  ) +
  theme_paper(base_size = 10.5)

save_plot(
  p_fig3,
  file.path(fig_main_dir, "03_fig3_null_calibration.pdf"),
  width = 8.6,
  height = 6.5
)

# ==============================================================================
# 5. Main-paper tables
# ==============================================================================

# ----------------------------------------------------------------------------
# Table 1: Detection algorithm summary
# ----------------------------------------------------------------------------

tab1_raw <- det_cell %>%
  group_by(method_id, method_label) %>%
  summarise(
    overlap = safe_mean(overlap_mean),
    overlap_q10 = safe_quantile(overlap_mean, 0.10),
    overlap_q90 = safe_quantile(overlap_mean, 0.90),
    ir = safe_median(influence_ratio_median),
    ir_q10 = safe_quantile(influence_ratio_median, 0.10),
    ir_q90 = safe_quantile(influence_ratio_median, 0.90),
    cpu = safe_median(cpu_median),
    cpu_q25 = safe_quantile(cpu_median, 0.25),
    cpu_q75 = safe_quantile(cpu_median, 0.75),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(method_id)

utils::write.csv(
  tab1_raw,
  file.path(data_dir, "03_tab1_algorithm_summary_numeric.csv"),
  row.names = FALSE,
  na = ""
)

tab1_display <- tab1_raw %>%
  transmute(
    `Detection method` = as.character(method_label),
    `Overlap, mean [P10, P90]` = fmt_pct_interval(
      overlap, overlap_q10, overlap_q90, digits = 1
    ),
    `Influence ratio, median [P10, P90]` = fmt_interval(
      ir, ir_q10, ir_q90, digits = 2
    ),
    `Runtime, median [IQR] (s)` = fmt_interval(
      cpu, cpu_q25, cpu_q75, digits = 4
    )
  )

write_tex_table(
  tab1_display,
  tex_path = file.path(tab_main_dir, "03_tab1_algorithm_summary.tex"),
  caption = paste0(
    "Detection accuracy, captured influence, and runtime across algorithmic ",
    "design cells. Brackets report the indicated between-cell percentiles."
  ),
  label = "tab:03-algorithm-summary",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrr"
)

# ----------------------------------------------------------------------------
# Table 2: EVT configuration summary under contaminated simulations
# ----------------------------------------------------------------------------

tab2_raw <- evt_cell %>%
  group_by(method_id, method_label) %>%
  summarise(
    convergence = safe_mean(convergence_rate),
    convergence_q10 = safe_quantile(convergence_rate, 0.10),
    convergence_q90 = safe_quantile(convergence_rate, 0.90),
    rejection = safe_mean(rejection_05),
    rejection_q10 = safe_quantile(rejection_05, 0.10),
    rejection_q90 = safe_quantile(rejection_05, 0.90),
    shape = safe_median(shape_median),
    shape_q25 = safe_quantile(shape_median, 0.25),
    shape_q75 = safe_quantile(shape_median, 0.75),
    cpu = safe_median(cpu_median),
    cpu_q25 = safe_quantile(cpu_median, 0.25),
    cpu_q75 = safe_quantile(cpu_median, 0.75),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(method_id)

utils::write.csv(
  tab2_raw,
  file.path(data_dir, "03_tab2_evt_summary_numeric.csv"),
  row.names = FALSE,
  na = ""
)

tab2_display <- tab2_raw %>%
  transmute(
    `EVT configuration` = as.character(method_label),
    `Convergence, mean [P10, P90]` = fmt_pct_interval(
      convergence, convergence_q10, convergence_q90, digits = 1
    ),
    `Rejection at 5%, mean [P10, P90]` = fmt_pct_interval(
      rejection, rejection_q10, rejection_q90, digits = 1
    ),
    `Shape, median [IQR]` = fmt_interval(
      shape, shape_q25, shape_q75, digits = 3
    ),
    `Runtime, median [IQR] (s)` = fmt_interval(
      cpu, cpu_q25, cpu_q75, digits = 4
    )
  )

write_tex_table(
  tab2_display,
  tex_path = file.path(tab_main_dir, "03_tab2_evt_configuration_summary.tex"),
  caption = paste0(
    "EVT convergence, signal rejection, fitted shape, and runtime under the ",
    "contaminated algorithmic simulations."
  ),
  label = "tab:03-evt-configuration-summary",
  resize_width = TABLE_WIDTH_WIDE,
  align = "lrrrr"
)

# ----------------------------------------------------------------------------
# Table 3: Null calibration summary
# ----------------------------------------------------------------------------

null_rate_summary <- null_alpha_cell %>%
  group_by(method_id, method_label, alpha) %>%
  summarise(
    empirical_size = safe_mean(rejection_rate),
    lower = safe_quantile(rejection_rate, 0.10),
    upper = safe_quantile(rejection_rate, 0.90),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # Use stable syntactic keys for pivot_wider(). Percentage-formatted
    # factor labels can vary across R/scales versions (for example "1%",
    # "1.0%", or strings containing padding), which previously caused the
    # expected "1%" column to be absent.
    alpha_key = case_when(
      dplyr::near(alpha, 0.010) ~ "alpha_001",
      dplyr::near(alpha, 0.025) ~ "alpha_0025",
      dplyr::near(alpha, 0.050) ~ "alpha_005",
      dplyr::near(alpha, 0.100) ~ "alpha_010",
      TRUE ~ NA_character_
    ),
    display = fmt_pct_interval(empirical_size, lower, upper, digits = 1)
  )

unknown_alpha_rows <- null_rate_summary %>%
  filter(is.na(alpha_key))

if (nrow(unknown_alpha_rows) > 0L) {
  warning(
    "Some null-calibration alpha values do not match ALPHA_GRID. ",
    "See diagnostics/03_unknown_null_alpha_values.csv."
  )
}

utils::write.csv(
  unknown_alpha_rows,
  file.path(diag_dir, "03_unknown_null_alpha_values.csv"),
  row.names = FALSE,
  na = ""
)

null_rate_summary <- null_rate_summary %>%
  filter(!is.na(alpha_key)) %>%
  select(method_id, method_label, alpha_key, display) %>%
  pivot_wider(
    names_from = alpha_key,
    values_from = display,
    names_expand = TRUE
  )

null_conv_summary <- null_cell %>%
  group_by(method_id, method_label) %>%
  summarise(
    convergence = safe_mean(convergence_rate),
    convergence_q10 = safe_quantile(convergence_rate, 0.10),
    convergence_q90 = safe_quantile(convergence_rate, 0.90),
    ks_fail_rate = safe_mean(ks_pvalue < 0.01),
    .groups = "drop"
  )

tab3_raw <- null_rate_summary %>%
  left_join(null_conv_summary, by = c("method_id", "method_label")) %>%
  arrange(method_id)

utils::write.csv(
  tab3_raw,
  file.path(data_dir, "03_tab3_null_calibration_numeric.csv"),
  row.names = FALSE,
  na = ""
)

# Ensure all expected alpha columns exist even when an entire method-alpha
# combination has no valid null-calibration observations.
expected_alpha_columns <- c(
  "alpha_001",
  "alpha_0025",
  "alpha_005",
  "alpha_010"
)

for (column_name in expected_alpha_columns) {
  if (!column_name %in% names(tab3_raw)) {
    tab3_raw[[column_name]] <- NA_character_
  }
}

tab3_display <- tab3_raw %>%
  transmute(
    `GEV fitting method` = as.character(method_label),
    `Nominal 1%` = .data[["alpha_001"]],
    `Nominal 2.5%` = .data[["alpha_0025"]],
    `Nominal 5%` = .data[["alpha_005"]],
    `Nominal 10%` = .data[["alpha_010"]],
    `Convergence, mean [P10, P90]` = fmt_pct_interval(
      convergence, convergence_q10, convergence_q90, digits = 1
    )
  )

write_tex_table(
  tab3_display,
  tex_path = file.path(tab_main_dir, "03_tab3_null_calibration.tex"),
  caption = paste0(
    "Empirical null rejection rates and convergence by GEV fitting method. ",
    "Rejection-rate brackets report the 10th and 90th percentiles across ",
    "architecture, sample-size, and set-size cells."
  ),
  label = "tab:03-null-calibration",
  resize_width = TABLE_WIDTH_WIDE,
  align = "lrrrrr"
)

# ==============================================================================
# 6. Supplementary figures
# ==============================================================================

# ----------------------------------------------------------------------------
# Figure A1: Detection overlap for every architecture
# ----------------------------------------------------------------------------

detection_architecture <- det_cell %>%
  group_by(architecture_label, k, method_id, method_label) %>%
  summarise(
    overlap = safe_mean(overlap_mean),
    lower = safe_quantile(overlap_mean, 0.10),
    upper = safe_quantile(overlap_mean, 0.90),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  filter(is.finite(overlap))

utils::write.csv(
  detection_architecture,
  file.path(data_dir, "03_figA1_detection_by_architecture_data.csv"),
  row.names = FALSE,
  na = ""
)

p_figA1 <- ggplot(
  detection_architecture,
  aes(
    x = k, y = overlap,
    colour = method_label, linetype = method_label,
    shape = method_label, group = method_label
  )
) +
  geom_line(linewidth = 0.62) +
  geom_point(size = 1.65, stroke = 0.65) +
  facet_wrap(~architecture_label, ncol = 3, scales = "free_x", drop = TRUE) +
  scale_colour_manual(values = det_method_colors, drop = FALSE) +
  scale_linetype_manual(values = det_method_linetypes, drop = FALSE) +
  scale_shape_manual(values = det_method_shapes, drop = FALSE) +
  scale_x_continuous(breaks = sort(unique(detection_architecture$k))) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  labs(
    x = "Influential-set size, k",
    y = "Mean detection overlap",
    colour = NULL, linetype = NULL, shape = NULL
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_paper(base_size = 9.4)

save_plot(
  p_figA1,
  file.path(fig_supp_dir, "03_figA1_detection_by_architecture.pdf"),
  width = 9.4,
  height = 10.8
)

# ----------------------------------------------------------------------------
# Figure A2: Exact-refined advantage over greedy-refined detection
# ----------------------------------------------------------------------------

algorithm_gain <- det_cell %>%
  filter(as.character(method_id) %in% c("greedy_ref", "exact_ref")) %>%
  select(
    architecture, architecture_label, N, k, B_type, B_actual, rho,
    method_id, overlap_mean
  ) %>%
  mutate(method_id = as.character(method_id)) %>%
  pivot_wider(names_from = method_id, values_from = overlap_mean) %>%
  mutate(gain = exact_ref - greedy_ref) %>%
  group_by(architecture, architecture_label, k) %>%
  summarise(
    gain = safe_mean(gain),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(k_label = factor(k, levels = sort(unique(k))))

max_gain <- max(abs(algorithm_gain$gain), na.rm = TRUE)
if (!is.finite(max_gain) || max_gain == 0) max_gain <- 0.01
algorithm_gain <- algorithm_gain %>%
  mutate(
    cell_label = sprintf("%+.1f", 100 * gain),
    text_color = ifelse(abs(gain) >= 0.58 * max_gain, "white", "black")
  )

utils::write.csv(
  algorithm_gain,
  file.path(data_dir, "03_figA2_exact_refined_advantage_data.csv"),
  row.names = FALSE,
  na = ""
)

p_figA2 <- ggplot(
  algorithm_gain,
  aes(x = k_label, y = architecture_label, fill = gain)
) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = cell_label, colour = text_color), size = 2.65) +
  scale_fill_gradient2(
    low = COL_ORANGE,
    mid = COL_NEUTRAL,
    high = COL_BLUE,
    midpoint = 0,
    limits = c(-max_gain, max_gain),
    oob = scales::squish,
    labels = scales::label_number(accuracy = 1, scale = 100, suffix = " pp"),
    name = "Exact-refined minus\ngreedy-refined overlap"
  ) +
  scale_colour_identity() +
  labs(
    x = "Influential-set size, k",
    y = NULL
  ) +
  guides(fill = guide_colourbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = grid::unit(7.5, "cm")
  )) +
  theme_heatmap(base_size = 9.8)

save_plot(
  p_figA2,
  file.path(fig_supp_dir, "03_figA2_exact_refined_advantage_heatmap.pdf"),
  width = 8.7,
  height = 6.4
)

# ----------------------------------------------------------------------------
# Figure A3: Calibration deviation at nominal 5%
# Orange = under-rejection; blue = over-rejection; white = nominal size.
# ----------------------------------------------------------------------------

calibration_heatmap <- null_alpha_cell %>%
  filter(abs(alpha - 0.05) < 1e-12) %>%
  group_by(architecture, architecture_label, k, method_id, method_label) %>%
  summarise(
    rejection_rate = safe_mean(rejection_rate),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(
    deviation = rejection_rate - 0.05,
    k_label = factor(k, levels = sort(unique(k)))
  )

max_cal_dev <- max(abs(calibration_heatmap$deviation), na.rm = TRUE)
if (!is.finite(max_cal_dev) || max_cal_dev == 0) max_cal_dev <- 0.01
calibration_heatmap <- calibration_heatmap %>%
  mutate(
    cell_label = sprintf("%.1f", 100 * rejection_rate),
    text_color = ifelse(abs(deviation) >= 0.58 * max_cal_dev, "white", "black")
  )

utils::write.csv(
  calibration_heatmap,
  file.path(data_dir, "03_figA3_calibration_heatmap_data.csv"),
  row.names = FALSE,
  na = ""
)

p_figA3 <- ggplot(
  calibration_heatmap,
  aes(x = k_label, y = architecture_label, fill = deviation)
) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = cell_label, colour = text_color), size = 2.65) +
  facet_wrap(~method_label, nrow = 1, drop = TRUE) +
  scale_fill_gradient2(
    low = COL_ORANGE,
    mid = COL_NEUTRAL,
    high = COL_BLUE,
    midpoint = 0,
    limits = c(-max_cal_dev, max_cal_dev),
    oob = scales::squish,
    labels = scales::label_number(accuracy = 1, scale = 100, suffix = " pp"),
    name = "Empirical size minus\nnominal 5%"
  ) +
  scale_colour_identity() +
  labs(
    x = "Set size, k",
    y = NULL
  ) +
  guides(fill = guide_colourbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = grid::unit(7.5, "cm")
  )) +
  theme_heatmap(base_size = 9.7)

save_plot(
  p_figA3,
  file.path(fig_supp_dir, "03_figA3_calibration_deviation_heatmap.pdf"),
  width = 9.3,
  height = 6.4
)

# ----------------------------------------------------------------------------
# Figure A4: Nestedness traces
# Open circles mark a violation of S_(k-1) being a subset of S_k.
# ----------------------------------------------------------------------------

if (nrow(nested) > 0L) {
  nested_method_labels <- c(
    "greedy" = "Greedy",
    "dinkelbach" = "Exact",
    "dinkelbach_refined" = "Exact + refinement"
  )
  nested_method_order <- c("greedy", "dinkelbach", "dinkelbach_refined")
  nested_colors <- c(
    "Greedy" = COL_GREY,
    "Exact" = COL_BLACK,
    "Exact + refinement" = COL_BLUE
  )
  nested_linetypes <- c(
    "Greedy" = "dotted",
    "Exact" = "dotdash",
    "Exact + refinement" = "solid"
  )
  nested_shapes <- c(
    "Greedy" = 1,
    "Exact" = 4,
    "Exact + refinement" = 16
  )
  
  nested_plot_data <- nested %>%
    mutate(
      method = factor(as.character(method), levels = nested_method_order),
      method_label = factor(
        unname(nested_method_labels[as.character(method)]),
        levels = unname(nested_method_labels[nested_method_order])
      ),
      violation = !is.na(nested) & !nested
    )
  
  utils::write.csv(
    nested_plot_data,
    file.path(data_dir, "03_figA4_nestedness_data.csv"),
    row.names = FALSE,
    na = ""
  )
  
  p_figA4 <- ggplot(
    nested_plot_data %>% filter(k > 1),
    aes(
      x = k, y = jaccard_with_prev,
      colour = method_label, linetype = method_label,
      shape = method_label, group = method_label
    )
  ) +
    geom_line(linewidth = 0.72) +
    geom_point(size = 2.0, stroke = 0.75) +
    geom_point(
      data = nested_plot_data %>% filter(k > 1, violation),
      shape = 21, fill = "white", size = 3.1, stroke = 0.9,
      show.legend = FALSE
    ) +
    facet_wrap(~architecture_label, ncol = 2, drop = TRUE) +
    scale_colour_manual(values = nested_colors, drop = FALSE) +
    scale_linetype_manual(values = nested_linetypes, drop = FALSE) +
    scale_shape_manual(values = nested_shapes, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0.01, 0.04))
    ) +
    labs(
      x = "Influential-set size, k",
      y = "Jaccard similarity with the preceding set",
      colour = NULL, linetype = NULL, shape = NULL
    ) +
    theme_paper(base_size = 10.2)
  
  save_plot(
    p_figA4,
    file.path(fig_supp_dir, "03_figA4_nestedness_traces.pdf"),
    width = 8.2,
    height = 6.3
  )
}

# ----------------------------------------------------------------------------
# Figure A5: Pooled null p-value QQ plots by architecture family
# ----------------------------------------------------------------------------

make_qq_data <- function(df) {
  p <- sort(df$p_value[is.finite(df$p_value)])
  n <- length(p)
  if (n == 0L) {
    return(data.frame(theoretical = numeric(0), empirical = numeric(0)))
  }
  
  # Retain at most 1,000 evenly spaced points per group for a compact PDF.
  index <- if (n <= 1000L) {
    seq_len(n)
  } else {
    unique(round(seq(1, n, length.out = 1000L)))
  }
  
  data.frame(
    theoretical = (index - 0.5) / n,
    empirical = p[index],
    stringsAsFactors = FALSE
  )
}

null_qq <- null_long %>%
  group_by(family_label, method_id, method_label) %>%
  group_modify(~make_qq_data(.x)) %>%
  ungroup()

utils::write.csv(
  null_qq,
  file.path(data_dir, "03_figA5_null_qq_data.csv"),
  row.names = FALSE,
  na = ""
)

p_figA5 <- ggplot(
  null_qq,
  aes(
    x = theoretical, y = empirical,
    colour = method_label, linetype = method_label,
    group = method_label
  )
) +
  geom_abline(
    intercept = 0, slope = 1,
    colour = COL_GREY, linewidth = 0.55, linetype = "longdash"
  ) +
  geom_line(linewidth = 0.75) +
  facet_wrap(~family_label, ncol = 2, drop = TRUE) +
  scale_colour_manual(values = null_method_colors, drop = FALSE) +
  scale_linetype_manual(values = null_method_linetypes, drop = FALSE) +
  scale_x_continuous(
    limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Uniform theoretical quantile",
    y = "Empirical p-value quantile",
    colour = NULL, linetype = NULL
  ) +
  coord_equal() +
  theme_paper(base_size = 10.1)

save_plot(
  p_figA5,
  file.path(fig_supp_dir, "03_figA5_null_pvalue_qq.pdf"),
  width = 8.0,
  height = 7.1
)

# ----------------------------------------------------------------------------
# Figure A6: Collinearity sensitivity
# Exact-refined minus greedy-refined overlap across rho, k, and N.
# ----------------------------------------------------------------------------

rho_gain <- det_cell %>%
  filter(
    architecture == "collinear_interaction",
    is.finite(rho),
    as.character(method_id) %in% c("greedy_ref", "exact_ref")
  ) %>%
  select(N, k, B_type, B_actual, rho, method_id, overlap_mean) %>%
  mutate(method_id = as.character(method_id)) %>%
  pivot_wider(names_from = method_id, values_from = overlap_mean) %>%
  mutate(gain = exact_ref - greedy_ref) %>%
  group_by(N, k, rho) %>%
  summarise(gain = safe_mean(gain), n_cells = n(), .groups = "drop") %>%
  mutate(
    N_label = factor(N, levels = sort(unique(N))),
    k_label = factor(k, levels = sort(unique(k))),
    rho_label = factor(rho, levels = sort(unique(rho)))
  )

if (nrow(rho_gain) > 0L) {
  max_rho_gain <- max(abs(rho_gain$gain), na.rm = TRUE)
  if (!is.finite(max_rho_gain) || max_rho_gain == 0) max_rho_gain <- 0.01
  
  rho_gain <- rho_gain %>%
    mutate(
      cell_label = sprintf("%+.1f", 100 * gain),
      text_color = ifelse(abs(gain) >= 0.58 * max_rho_gain, "white", "black")
    )
  
  utils::write.csv(
    rho_gain,
    file.path(data_dir, "03_figA6_collinearity_sensitivity_data.csv"),
    row.names = FALSE,
    na = ""
  )
  
  p_figA6 <- ggplot(
    rho_gain,
    aes(x = rho_label, y = k_label, fill = gain)
  ) +
    geom_tile(colour = "white", linewidth = 0.45) +
    geom_text(aes(label = cell_label, colour = text_color), size = 2.55) +
    facet_wrap(~N_label, nrow = 1, labeller = label_both) +
    scale_fill_gradient2(
      low = COL_ORANGE,
      mid = COL_NEUTRAL,
      high = COL_BLUE,
      midpoint = 0,
      limits = c(-max_rho_gain, max_rho_gain),
      oob = scales::squish,
      labels = scales::label_number(accuracy = 1, scale = 100, suffix = " pp"),
      name = "Exact-refined minus\ngreedy-refined overlap"
    ) +
    scale_colour_identity() +
    labs(
      x = "Predictor correlation, rho",
      y = "Influential-set size, k"
    ) +
    guides(fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(7.5, "cm")
    )) +
    theme_heatmap(base_size = 9.6)
  
  save_plot(
    p_figA6,
    file.path(fig_supp_dir, "03_figA6_collinearity_sensitivity.pdf"),
    width = 9.4,
    height = 4.8
  )
}

# ----------------------------------------------------------------------------
# Figure A7: EVT convergence by architecture and sample size
# Sequential blue is used because convergence is not a signed contrast.
# ----------------------------------------------------------------------------

evt_convergence_heatmap <- evt_cell %>%
  group_by(architecture, architecture_label, N, method_id, method_label) %>%
  summarise(
    convergence_rate = safe_mean(convergence_rate),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(
    N_label = factor(N, levels = sort(unique(N))),
    cell_label = sprintf("%.0f", 100 * convergence_rate),
    text_color = ifelse(convergence_rate >= 0.62, "white", "black")
  )

utils::write.csv(
  evt_convergence_heatmap,
  file.path(data_dir, "03_figA7_evt_convergence_data.csv"),
  row.names = FALSE,
  na = ""
)

p_figA7 <- ggplot(
  evt_convergence_heatmap,
  aes(x = N_label, y = architecture_label, fill = convergence_rate)
) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = cell_label, colour = text_color), size = 2.55) +
  facet_wrap(~method_label, nrow = 1, drop = TRUE) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = COL_BLUE_DARK,
    limits = c(0, 1),
    oob = scales::squish,
    labels = scales::label_percent(accuracy = 1),
    name = "Convergence rate"
  ) +
  scale_colour_identity() +
  labs(
    x = "Sample size",
    y = NULL
  ) +
  guides(fill = guide_colourbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = grid::unit(7.0, "cm")
  )) +
  theme_heatmap(base_size = 9.3)

save_plot(
  p_figA7,
  file.path(fig_supp_dir, "03_figA7_evt_convergence_heatmap.pdf"),
  width = 10.1,
  height = 6.3
)

# ==============================================================================
# 7. Supplementary tables
# ==============================================================================

# ----------------------------------------------------------------------------
# Table A1: Architecture-level refined-method comparison
# ----------------------------------------------------------------------------

tabA1_raw <- det_cell %>%
  filter(as.character(method_id) %in% c("greedy_ref", "exact_ref")) %>%
  group_by(architecture, architecture_label, method_id) %>%
  summarise(
    overlap = safe_mean(overlap_mean),
    ir = safe_median(influence_ratio_median),
    cpu = safe_median(cpu_median),
    .groups = "drop"
  ) %>%
  mutate(method_id = as.character(method_id)) %>%
  pivot_wider(
    names_from = method_id,
    values_from = c(overlap, ir, cpu)
  ) %>%
  mutate(
    overlap_gain = overlap_exact_ref - overlap_greedy_ref,
    runtime_ratio = cpu_exact_ref / cpu_greedy_ref
  ) %>%
  arrange(architecture_label)

utils::write.csv(
  tabA1_raw,
  file.path(data_dir, "03_tabA1_architecture_comparison_numeric.csv"),
  row.names = FALSE,
  na = ""
)

tabA1_display <- tabA1_raw %>%
  transmute(
    Architecture = as.character(architecture_label),
    `Greedy-refined overlap` = fmt_pct(overlap_greedy_ref),
    `Exact-refined overlap` = fmt_pct(overlap_exact_ref),
    `Exact advantage` = fmt_pp(overlap_gain),
    `Greedy-refined influence ratio` = fmt_num(ir_greedy_ref, 2),
    `Exact-refined influence ratio` = fmt_num(ir_exact_ref, 2),
    `Exact/greedy runtime` = fmt_num(runtime_ratio, 2)
  )

write_tex_table(
  tabA1_display,
  tex_path = file.path(tab_supp_dir, "03_tabA1_architecture_comparison.tex"),
  caption = paste0(
    "Architecture-level comparison of the refined greedy and refined exact ",
    "detection methods."
  ),
  label = "tab:03A-architecture-comparison",
  resize_width = TABLE_WIDTH_WIDE,
  align = "lrrrrrr",
  font_command = "\\small"
)

# ----------------------------------------------------------------------------
# Table A2: Null calibration at 5% by architecture
# ----------------------------------------------------------------------------

tabA2_rates <- null_alpha_cell %>%
  filter(abs(alpha - 0.05) < 1e-12) %>%
  group_by(architecture, architecture_label, method_id, method_label) %>%
  summarise(
    size = safe_mean(rejection_rate),
    size_q10 = safe_quantile(rejection_rate, 0.10),
    size_q90 = safe_quantile(rejection_rate, 0.90),
    .groups = "drop"
  )

tabA2_conv <- null_cell %>%
  group_by(architecture, architecture_label, method_id, method_label) %>%
  summarise(
    convergence = safe_mean(convergence_rate),
    ks_fail_rate = safe_mean(ks_pvalue < 0.01),
    shape = safe_median(shape_median),
    .groups = "drop"
  )

tabA2_raw <- tabA2_rates %>%
  left_join(
    tabA2_conv,
    by = c("architecture", "architecture_label", "method_id", "method_label")
  ) %>%
  arrange(architecture_label, method_id)

utils::write.csv(
  tabA2_raw,
  file.path(data_dir, "03_tabA2_calibration_by_architecture_numeric.csv"),
  row.names = FALSE,
  na = ""
)

tabA2_display <- tabA2_raw %>%
  transmute(
    Architecture = as.character(architecture_label),
    `GEV fit` = as.character(method_label),
    `Size at 5%, mean [P10, P90]` = fmt_pct_interval(
      size, size_q10, size_q90, digits = 1
    ),
    Convergence = fmt_pct(convergence),
    `KS p below 0.01` = fmt_pct(ks_fail_rate),
    `Median shape` = fmt_num(shape, 3)
  )

write_tex_table(
  tabA2_display,
  tex_path = file.path(tab_supp_dir, "03_tabA2_calibration_by_architecture.tex"),
  caption = paste0(
    "Null calibration, convergence, and fitted shape by architecture and GEV ",
    "fitting method."
  ),
  label = "tab:03A-calibration-by-architecture",
  resize_width = TABLE_WIDTH_WIDE,
  align = "llrrrr",
  font_command = "\\small"
)

# ----------------------------------------------------------------------------
# Table A3: Adaptive block-count capping
# Count each scenario once rather than once per Monte Carlo iteration.
# ----------------------------------------------------------------------------

tabA3_raw <- scaling %>%
  distinct(
    architecture, architecture_label, N, k, B_type,
    B_raw, B_actual, B_capped, rho
  ) %>%
  group_by(N, B_type) %>%
  summarise(
    scenarios = n(),
    capping_rate = safe_mean(as.numeric(B_capped)),
    median_requested = safe_median(B_raw),
    median_actual = safe_median(B_actual),
    minimum_actual = min(B_actual, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(N, B_type)

utils::write.csv(
  tabA3_raw,
  file.path(data_dir, "03_tabA3_block_capping_numeric.csv"),
  row.names = FALSE,
  na = ""
)

tabA3_display <- tabA3_raw %>%
  transmute(
    `Sample size` = N,
    `Requested block rule` = B_type,
    Scenarios = scenarios,
    `Capping rate` = fmt_pct(capping_rate),
    `Median requested` = fmt_num(median_requested, 0),
    `Median actual` = fmt_num(median_actual, 0),
    `Minimum actual` = minimum_actual
  )

write_tex_table(
  tabA3_display,
  tex_path = file.path(tab_supp_dir, "03_tabA3_block_capping.tex"),
  caption = "Adaptive block-count capping by sample size and requested block rule.",
  label = "tab:03A-block-capping",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "llrrrrr"
)

# ----------------------------------------------------------------------------
# Table A4: EVT performance by architecture
# ----------------------------------------------------------------------------

tabA4_raw <- evt_cell %>%
  group_by(architecture, architecture_label, method_id, method_label) %>%
  summarise(
    convergence = safe_mean(convergence_rate),
    rejection = safe_mean(rejection_05),
    shape = safe_median(shape_median),
    cpu = safe_median(cpu_median),
    .groups = "drop"
  ) %>%
  arrange(architecture_label, method_id)

utils::write.csv(
  tabA4_raw,
  file.path(data_dir, "03_tabA4_evt_by_architecture_numeric.csv"),
  row.names = FALSE,
  na = ""
)

tabA4_display <- tabA4_raw %>%
  transmute(
    Architecture = as.character(architecture_label),
    `EVT configuration` = as.character(method_label),
    Convergence = fmt_pct(convergence),
    `Rejection at 5%` = fmt_pct(rejection),
    `Median shape` = fmt_num(shape, 3),
    `Median runtime (s)` = fmt_num(cpu, 4)
  )

write_tex_table(
  tabA4_display,
  tex_path = file.path(tab_supp_dir, "03_tabA4_evt_by_architecture.tex"),
  caption = "EVT performance by architecture and fitting configuration.",
  label = "tab:03A-evt-by-architecture",
  resize_width = TABLE_WIDTH_WIDE,
  align = "llrrrr",
  font_command = "\\small"
)

# ----------------------------------------------------------------------------
# Table A5: Nestedness summary
# ----------------------------------------------------------------------------

if (nrow(nested) > 0L) {
  nested_method_labels <- c(
    "greedy" = "Greedy",
    "dinkelbach" = "Exact",
    "dinkelbach_refined" = "Exact + refinement"
  )
  
  tabA5_raw <- nested %>%
    mutate(
      method_label = unname(nested_method_labels[as.character(method)])
    ) %>%
    group_by(architecture, architecture_label, method, method_label) %>%
    summarise(
      mean_jaccard = safe_mean(jaccard_with_prev[k > 1]),
      minimum_jaccard = {
        x <- jaccard_with_prev[k > 1 & is.finite(jaccard_with_prev)]
        if (length(x) == 0L) NA_real_ else min(x)
      },
      violations = sum(nested %in% FALSE, na.rm = TRUE),
      median_cpu = safe_median(cpu_seconds[cpu_seconds > 0]),
      final_influence = {
        valid <- which(is.finite(influence_magnitude))
        if (length(valid) == 0L) NA_real_ else influence_magnitude[max(valid)]
      },
      .groups = "drop"
    ) %>%
    arrange(architecture_label, method)
  
  utils::write.csv(
    tabA5_raw,
    file.path(data_dir, "03_tabA5_nestedness_summary_numeric.csv"),
    row.names = FALSE,
    na = ""
  )
  
  tabA5_display <- tabA5_raw %>%
    transmute(
      Architecture = as.character(architecture_label),
      Method = method_label,
      `Mean Jaccard` = fmt_num(mean_jaccard, 3),
      `Minimum Jaccard` = fmt_num(minimum_jaccard, 3),
      Violations = violations,
      `Median time (s)` = fmt_num(median_cpu, 4),
      `Final influence` = fmt_num(final_influence, 3)
    )
  
  write_tex_table(
    tabA5_display,
    tex_path = file.path(tab_supp_dir, "03_tabA5_nestedness_summary.tex"),
    caption = "Nestedness and computational summaries for the representative traces.",
    label = "tab:03A-nestedness-summary",
    resize_width = TABLE_WIDTH_WIDE,
    align = "llrrrrr"
  )
}

# Full design-cell results are retained as CSV because complete LaTeX versions
# would be too tall for a formal single-page table environment.
utils::write.csv(
  det_cell,
  file.path(tab_supp_dir, "03_tabA6_full_detection_design_cells.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  evt_cell,
  file.path(tab_supp_dir, "03_tabA7_full_evt_design_cells.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  null_alpha_cell,
  file.path(tab_supp_dir, "03_tabA8_full_null_calibration_cells.csv"),
  row.names = FALSE,
  na = ""
)

# ==============================================================================
# 8. Diagnostics and manifest
# ==============================================================================

# Input and failure audit.
input_audit <- data.frame(
  dataset = c(
    "03v2 scaling", "03b null calibration",
    if (nrow(nested) > 0L) "03v2 nestedness" else character(0)
  ),
  rows = c(
    nrow(scaling), nrow(null_cal),
    if (nrow(nested) > 0L) nrow(nested) else integer(0)
  ),
  columns = c(
    ncol(scaling), ncol(null_cal),
    if (nrow(nested) > 0L) ncol(nested) else integer(0)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  input_audit,
  file.path(diag_dir, "03_input_audit.csv"),
  row.names = FALSE
)

scaling_failures <- scaling %>%
  mutate(error_msg = as.character(error_msg)) %>%
  filter(!is.na(error_msg), nzchar(error_msg)) %>%
  count(error_msg, sort = TRUE, name = "n")
utils::write.csv(
  scaling_failures,
  file.path(diag_dir, "03_scaling_failure_messages.csv"),
  row.names = FALSE,
  na = ""
)

null_failures <- null_cal %>%
  mutate(error_msg = as.character(error_msg)) %>%
  filter(!is.na(error_msg), nzchar(error_msg)) %>%
  count(error_msg, sort = TRUE, name = "n")
utils::write.csv(
  null_failures,
  file.path(diag_dir, "03_null_failure_messages.csv"),
  row.names = FALSE,
  na = ""
)

capture.output(
  sessionInfo(),
  file = file.path(diag_dir, "03_output_session_info.txt")
)

manifest_path <- file.path(diag_dir, "03_output_manifest.csv")
all_outputs <- list.files(
  output_root,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)
all_outputs <- setdiff(
  normalizePath(all_outputs, winslash = "/", mustWork = FALSE),
  normalizePath(manifest_path, winslash = "/", mustWork = FALSE)
)
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

cat("\nScript 83 completed successfully.\n")
cat("Output directory:\n  ", output_root, "\n", sep = "")
cat(sprintf("Generated %d output files.\n", nrow(manifest) + 1L))
