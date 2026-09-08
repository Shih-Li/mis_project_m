# ==============================================================================
# File: scripts/84_output_robust.R
# Purpose:
#   Generate publication-ready figures and tables for Script 04:
#   comparison of OLS, classical diagnostic deletion, MIS with oracle k,
#   MM regression, and LTS, with emphasis on influential-set recovery
#   under bad-leverage contamination.
#
# Primary input:
#   output/04sap_robust_comparison_results.rds
#
# Optional validation inputs:
#   output/04sap_summary_tables.rds
#   output/04sap_bias_distributional.rds
#
# Output structure:
#   output/04_robust/
#     figures/main/          Main-paper PDF figures
#     figures/supplement/    Supplementary PDF figures
#     tables/main/           Complete LaTeX table environments + CSV companions
#     tables/supplement/     Supplementary LaTeX tables + CSV companions
#     data/                  Data used to construct figures and tables
#     diagnostics/           Audits, session information, and output manifest
#
# LaTeX requirements for generated tables:
#   \usepackage{graphicx}
#   \usepackage{booktabs}
#
# Statistical reporting rules:
#   - Arithmetic means are the primary summaries.
#   - No median-based headline results are produced.
#   - Mean coefficient, mean absolute bias, RMSE, empirical coverage,
#     mean selected k, mean overlap, and mean runtime are reported.
#   - Monte Carlo standard errors are calculated whenever appropriate.
#   - Main broad comparisons give equal weight to each simulation design cell.
#
# Figure design rules:
#   - No title or subtitle inside plots; captions are handled in LaTeX.
#   - Full empirical coefficient-error distributions are shown using
#   iteration-level points on a signed pseudo-logarithmic scale.
#   - Arithmetic means are the primary summaries.
#   - Medians are shown only as secondary markers of skewness and
#   are not used as headline performance measures.
#   - Labels, legends, facets, and margins are sized for formal papers.
#   - No red-green palette is used.
#   - Large-error heatmaps use blue-white-orange:
#          blue = lower exceedance probability;
#          orange = higher exceedance probability.
#   - Main Figure 2 focuses on MIS with oracle k.
#   - The all-estimator large-error heatmap is supplementary.
#
# Table design rules:
#   - Tables are complete \begin{table}[htbp] floating environments.
#   - Compact tables use 0.7\columnwidth.
#   - Wider tables use 0.85\columnwidth or \columnwidth when necessary
#     to preserve readable labels and indicators.
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
    ". Install them before running scripts/84_output_robust.R."
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

input_main <- file.path(
  project_root, "output", "04sap_robust_comparison_results.rds"
)
input_summary_optional <- file.path(
  project_root, "output", "04sap_summary_tables.rds"
)
input_bias_optional <- file.path(
  project_root, "output", "04sap_bias_distributional.rds"
)

output_root <- file.path(project_root, "output", "04_robust")

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

if (!file.exists(input_main)) {
  stop("Missing primary input file: ", input_main)
}


# Publication output.
SAVE_PNG_PREVIEWS <- FALSE
PNG_DPI <- 320

# True coefficient used by the current Script 04 simulation.
TRUE_BETA <- 1

# Color-blind-conscious palette: no red-green scale.
COL_ORANGE      <- "#E69F00"
COL_ORANGE_DARK <- "#B35806"
COL_BLUE        <- "#0072B2"
COL_BLUE_LIGHT  <- "#56B4E9"
COL_BLUE_DARK   <- "#08519C"
COL_PURPLE      <- "#6A3D9A"
COL_BLACK       <- "#111111"
COL_GREY_DARK   <- "#4D4D4D"
COL_GREY        <- "#7A7A7A"
COL_GREY_LIGHT  <- "#B8B8B8"
COL_NEUTRAL     <- "#F7F7F7"
COL_VIOLIN      <- "#DCE6EF"

TABLE_WIDTH_COMPACT <- "0.7\\columnwidth"
TABLE_WIDTH_MEDIUM  <- "0.85\\columnwidth"
TABLE_WIDTH_WIDE    <- "\\columnwidth"


# ==============================================================================
# 1. Shared helper functions
# ==============================================================================

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}


safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}


safe_mcse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}


safe_prop_mcse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  p_hat <- mean(x)
  sqrt(p_hat * (1 - p_hat) / length(x))
}


safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  unname(stats::quantile(
    x,
    probs = probability,
    na.rm = TRUE,
    names = FALSE,
    type = 7
  ))
}


combined_cell_mcse <- function(cell_mcse, cell_mean) {
  valid <- is.finite(cell_mcse) & is.finite(cell_mean)
  n_cells <- sum(valid)
  
  if (n_cells == 0L) return(NA_real_)
  
  sqrt(sum(cell_mcse[valid]^2)) / n_cells
}


fmt_num <- function(x, digits = 3L) {
  ifelse(
    is.na(x),
    "--",
    formatC(x, format = "f", digits = digits)
  )
}


fmt_pct <- function(x, digits = 1L) {
  ifelse(
    is.na(x),
    "--",
    paste0(formatC(100 * x, format = "f", digits = digits), "%")
  )
}


fmt_mean_mcse <- function(mean_value, mcse_value, digits = 3L) {
  ifelse(
    is.na(mean_value),
    "--",
    paste0(
      formatC(mean_value, format = "f", digits = digits),
      " (",
      ifelse(
        is.na(mcse_value),
        "--",
        formatC(mcse_value, format = "f", digits = digits)
      ),
      ")"
    )
  )
}


fmt_pct_mcse <- function(mean_value, mcse_value, digits = 1L) {
  ifelse(
    is.na(mean_value),
    "--",
    paste0(
      formatC(100 * mean_value, format = "f", digits = digits),
      "\\% (",
      ifelse(
        is.na(mcse_value),
        "--",
        formatC(100 * mcse_value, format = "f", digits = digits)
      ),
      ")"
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
    placement = "htbp",
    font_command = NULL,
    escape_cells = TRUE
) {
  if (!is.data.frame(data) || ncol(data) == 0L) {
    stop("write_tex_table() requires a non-empty data.frame.")
  }
  
  utils::write.csv(data, csv_path, row.names = FALSE, na = "")
  
  if (is.null(align)) {
    align <- paste0(
      "l",
      paste(rep("r", max(0L, ncol(data) - 1L)), collapse = "")
    )
  }
  
  if (nchar(align) != ncol(data)) {
    stop("Length of LaTeX alignment string must equal the number of columns.")
  }
  
  escaped_names <- escape_latex(names(data))
  
  if (isTRUE(escape_cells)) {
    escaped_data <- lapply(data, escape_latex)
  } else {
    escaped_data <- lapply(data, as.character)
  }
  
  escaped_data <- as.data.frame(
    escaped_data,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
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
      plot.margin = margin(8, 12, 8, 8)
    )
}


theme_distribution <- function(base_size = 10.5, base_family = "sans") {
  theme_paper(base_size = base_size, base_family = base_family) +
    theme(
      axis.text.y = element_text(size = base_size - 0.2),
      panel.grid.major.x = element_line(
        colour = "grey90",
        linewidth = 0.25
      ),
      panel.grid.minor = element_blank()
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
      plot.margin = margin(8, 12, 8, 8)
    )
}


reshape_mapped_columns <- function(
    data,
    mapping,
    value_name,
    metadata,
    metadata_key
) {
  scenario_columns <- c(
    "iter",
    "n_obs",
    "design_k",
    "contam_prop",
    "x_type",
    "error_type",
    "outlier_method",
    "set_size"
  )
  
  result <- data %>%
    select(all_of(scenario_columns), all_of(unname(mapping))) %>%
    pivot_longer(
      cols = all_of(unname(mapping)),
      names_to = "source_column",
      values_to = value_name
    ) %>%
    mutate(
      !!metadata_key := names(mapping)[
        match(source_column, unname(mapping))
      ]
    ) %>%
    select(-source_column)
  
  left_join(
    result,
    metadata,
    by = setNames(metadata_key, metadata_key)
  )
}


add_outlier_display <- function(data) {
  data %>%
    mutate(
      outlier_label_plot = factor(
        unname(outlier_labels_plot[as.character(outlier_method)]),
        levels = unname(outlier_labels_plot[outlier_order])
      ),
      outlier_label_table = factor(
        unname(outlier_labels_table[as.character(outlier_method)]),
        levels = unname(outlier_labels_table[outlier_order])
      )
    )
}


make_metric_wide_table <- function(
    data,
    mean_column,
    mcse_column,
    formatter,
    include_clean = TRUE
) {
  working <- data
  
  if (!isTRUE(include_clean)) {
    working <- working %>%
      filter(outlier_method != "none")
  }
  
  working %>%
    add_outlier_display() %>%
    mutate(
      cell = formatter(
        .data[[mean_column]],
        .data[[mcse_column]]
      )
    ) %>%
    select(
      estimator_order,
      Estimator = estimator_label,
      outlier_label_table,
      cell
    ) %>%
    arrange(estimator_order, outlier_label_table) %>%
    select(-estimator_order) %>%
    pivot_wider(
      names_from = outlier_label_table,
      values_from = cell
    ) %>%
    select(
      Estimator,
      any_of(unname(
        outlier_labels_table[
          if (include_clean) outlier_order else outlier_order[-1L]
        ]
      ))
    )
}


# ==============================================================================
# 2. Load and validate inputs
# ==============================================================================

sim <- readRDS(input_main)


# ==============================================================================
# Legacy Script 04 compatibility
# ==============================================================================

# The legacy Script 04 simulation used n = 5000 and stored set_size,
# but did not store n_obs, design_k, or contam_prop explicitly.

if (!"n_obs" %in% names(sim)) {
  sim$n_obs <- 5000L
}

if (!"design_k" %in% names(sim)) {
  sim$design_k <- as.integer(sim$set_size)
}

if (!"contam_prop" %in% names(sim)) {
  sim$contam_prop <- ifelse(
    sim$outlier_method == "none",
    0,
    sim$design_k / sim$n_obs
  )
}


# ==============================================================================
# GPD handling
# ==============================================================================

# Retain all Script 04 conditions for condition-specific diagnostics.
sim_all <- sim

# Primary finite-moment comparison excludes the legacy GPD specification.
sim_primary <- sim %>%
  filter(error_type != "gpd")

n_gpd <- sum(
  sim_all$error_type == "gpd",
  na.rm = TRUE
)

cat(sprintf(
  "GPD rows retained for condition-specific diagnostics: %s\n",
  format(n_gpd, big.mark = ",")
))

required_columns <- c(
  "iter",
  "n_obs",
  "design_k",
  "contam_prop",
  "x_type",
  "error_type",
  "outlier_method",
  "set_size",
  
  "k_cd", "k_lev", "k_dfb", "k_alpha", "k_oracle",
  "k_peel_v2", "k_peel_sap",
  
  "overlap_cd", "overlap_lev", "overlap_dfb",
  "overlap_mis_alpha", "overlap_mis_oracle",
  "overlap_peel_v2", "overlap_peel_sap",
  
  "peel_v2_stop", "peel_v2_iters",
  "peel_sap_stop", "peel_sap_iters",
  "peel_sap_final_p", "peel_sap_min_p",
  "peel_sap_direction", "peel_sap_peak_excess",
  
  "mm_converged", "mm_valid", "mm_zero_scale", "mm_scale",
  
  "coef_full", "coef_cd", "coef_lev", "coef_dfb",
  "coef_mis_alpha", "coef_mis_oracle",
  "coef_mis_peel", "coef_mis_sap",
  "coef_mm", "coef_lts",
  
  "se_full", "se_cd", "se_lev", "se_dfb",
  "se_mis_alpha", "se_mis_oracle",
  "se_mis_peel", "se_mis_sap",
  "se_mm", "se_lts",
  
  "bias_full", "bias_cd", "bias_lev", "bias_dfb",
  "bias_mis_alpha", "bias_mis_oracle",
  "bias_mis_peel", "bias_mis_sap",
  "bias_mm", "bias_lts",
  
  "cov_full", "cov_cd", "cov_lev", "cov_dfb",
  "cov_mis_alpha", "cov_mis_oracle",
  "cov_mis_peel", "cov_mis_sap",
  "cov_mm", "cov_lts",
  
  "cpu_full", "cpu_cd", "cpu_lev", "cpu_dfb",
  "cpu_mis_alpha", "cpu_mis_oracle",
  "cpu_peel_v2", "cpu_peel_sap",
  "cpu_mm", "cpu_lts"
)

missing_columns <- setdiff(required_columns, names(sim))

if (length(missing_columns) > 0L) {
  stop(
    "04sap_robust_comparison_results.rds is missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

n_obs_levels <- sim %>%
  distinct(n_obs) %>%
  arrange(n_obs) %>%
  pull(n_obs)

contam_prop_levels <- sim %>%
  filter(outlier_method != "none") %>%
  distinct(contam_prop) %>%
  arrange(contam_prop) %>%
  pull(contam_prop)

expected_n_obs_levels <- c(
  500L,
  1000L,
  2500L,
  5000L
)

expected_contam_prop_levels <- c(
  0.005,
  0.010,
  0.025,
  0.050
)

if (!all(
  as.integer(n_obs_levels) %in%
  expected_n_obs_levels
)) {
  warning(
    "Unexpected sample-size grid. Found: ",
    paste(n_obs_levels, collapse = ", ")
  )
}

if (!all(
  as.numeric(contam_prop_levels) %in%
  expected_contam_prop_levels
)) {
  warning(
    "Unexpected contamination-proportion grid. Found: ",
    paste(contam_prop_levels, collapse = ", ")
  )
}

cat(sprintf(
  "Loaded Script 04 results: %s rows\n",
  format(nrow(sim), big.mark = ",")
))

cat(sprintf(
  "Optional summary RDS present: %s\n",
  ifelse(file.exists(input_summary_optional), "yes", "no")
))

cat(sprintf(
  "Optional bias-summary RDS present: %s\n",
  ifelse(file.exists(input_bias_optional), "yes", "no")
))


# ==============================================================================
# 3. Labels, ordering, and method mappings
# ==============================================================================

outlier_order <- c(
  "none",
  "vertical_outlier",
  "good_leverage",
  "bad_leverage"
)

outlier_labels_plot <- c(
  "none" = "No contamination",
  "vertical_outlier" = "Vertical outliers",
  "good_leverage" = "Good leverage",
  "bad_leverage" = "Bad leverage"
)

outlier_labels_table <- c(
  "none" = "None",
  "vertical_outlier" = "Vertical",
  "good_leverage" = "Good leverage",
  "bad_leverage" = "Bad leverage"
)

error_order <- c(
  "normal",
  "mixed_normal",
  "beta_logistic",
  "skewed_t",
  "contaminated",
  "golm",
  "pareto",
  "gpd"
)

error_labels_plot <- c(
  "normal" = "Normal",
  "mixed_normal" = "Mixed\nnormal",
  "beta_logistic" = "Beta-\nlogistic",
  "skewed_t" = "Skewed-t",
  "contaminated" = "Contaminated",
  "golm" = "GOLM",
  "pareto" = "Pareto",
  "gpd" = "GPD"
)

error_labels_table <- c(
  "normal" = "Normal",
  "mixed_normal" = "Mixed normal",
  "beta_logistic" = "Beta-logistic",
  "skewed_t" = "Skewed-t",
  "contaminated" = "Contaminated",
  "golm" = "GOLM",
  "pareto" = "Pareto",
  "gpd" = "GPD"
)

x_order <- c("normal", "mixed_normal", "contaminated")

x_labels_plot <- c(
  "normal" = "Normal",
  "mixed_normal" = "Mixed normal",
  "contaminated" = "Contaminated"
)

x_labels_table <- c(
  "normal" = "Normal",
  "mixed_normal" = "Mixed normal",
  "contaminated" = "Contaminated"
)


estimator_meta <- data.frame(
  estimator_id = c(
    "full",
    "cd",
    "lev",
    "dfb",
    "mis_alpha",
    "mis_peel",
    "mis_sap",
    "mis_oracle",
    "mm",
    "lts"
  ),
  estimator_label = c(
    "OLS",
    "Cook's D",
    "Leverage",
    "DFBETAS",
    "MIS, adaptive k",
    "MIS, iterative peel",
    "MIS-SAP",
    "MIS, oracle k",
    "MM-estimator",
    "LTS"
  ),
  estimator_short = c(
    "OLS",
    "CD",
    "LEV",
    "DFB",
    "MIS-a",
    "Peel",
    "SAP",
    "Oracle",
    "MM",
    "LTS"
  ),
  estimator_family = c(
    "OLS",
    "Classical deletion",
    "Classical deletion",
    "Classical deletion",
    "Data-driven MIS",
    "Data-driven MIS",
    "Data-driven MIS",
    "Oracle benchmark",
    "Robust regression",
    "Robust regression"
  ),
  estimator_order = seq_len(10L),
  stringsAsFactors = FALSE
)


selection_meta <- estimator_meta %>%
  filter(estimator_id %in% c(
    "cd", "lev", "dfb",
    "mis_alpha", "mis_peel", "mis_sap", "mis_oracle"
  )) %>%
  mutate(
    selection_order = match(
      estimator_id,
      c(
        "cd", "lev", "dfb",
        "mis_alpha", "mis_peel", "mis_sap", "mis_oracle"
      )
    )
  )


method_colors <- c(
  "full" = COL_BLACK,
  "cd" = "#969696",
  "lev" = "#BDBDBD",
  "dfb" = "#636363",
  "mis_alpha" = COL_BLUE_LIGHT,
  "mis_peel" = COL_BLUE,
  "mis_sap" = COL_BLUE_DARK,
  "mis_oracle" = COL_PURPLE,
  "mm" = COL_ORANGE,
  "lts" = COL_ORANGE_DARK
)

method_shapes <- c(
  "full" = 16,
  "cd" = 0,
  "lev" = 1,
  "dfb" = 2,
  "mis_alpha" = 15,
  "mis_peel" = 17,
  "mis_sap" = 18,
  "mis_oracle" = 8,
  "mm" = 3,
  "lts" = 4
)

method_linetypes <- c(
  "full" = "solid",
  "cd" = "dotted",
  "lev" = "dotdash",
  "dfb" = "longdash",
  "mis_alpha" = "dashed",
  "mis_peel" = "solid",
  "mis_sap" = "twodash",
  "mis_oracle" = "longdash",
  "mm" = "dashed",
  "lts" = "dotdash"
)


coef_mapping <- c(
  "full" = "coef_full",
  "cd" = "coef_cd",
  "lev" = "coef_lev",
  "dfb" = "coef_dfb",
  "mis_alpha" = "coef_mis_alpha",
  "mis_peel" = "coef_mis_peel",
  "mis_sap" = "coef_mis_sap",
  "mis_oracle" = "coef_mis_oracle",
  "mm" = "coef_mm",
  "lts" = "coef_lts"
)

bias_mapping <- c(
  "full" = "bias_full",
  "cd" = "bias_cd",
  "lev" = "bias_lev",
  "dfb" = "bias_dfb",
  "mis_alpha" = "bias_mis_alpha",
  "mis_peel" = "bias_mis_peel",
  "mis_sap" = "bias_mis_sap",
  "mis_oracle" = "bias_mis_oracle",
  "mm" = "bias_mm",
  "lts" = "bias_lts"
)

coverage_mapping <- c(
  "full" = "cov_full",
  "cd" = "cov_cd",
  "lev" = "cov_lev",
  "dfb" = "cov_dfb",
  "mis_alpha" = "cov_mis_alpha",
  "mis_peel" = "cov_mis_peel",
  "mis_sap" = "cov_mis_sap",
  "mis_oracle" = "cov_mis_oracle",
  "mm" = "cov_mm",
  "lts" = "cov_lts"
)

runtime_mapping <- c(
  "full" = "cpu_full",
  "cd" = "cpu_cd",
  "lev" = "cpu_lev",
  "dfb" = "cpu_dfb",
  "mis_alpha" = "cpu_mis_alpha",
  "mis_peel" = "cpu_peel_v2",
  "mis_sap" = "cpu_peel_sap",
  "mis_oracle" = "cpu_mis_oracle",
  "mm" = "cpu_mm",
  "lts" = "cpu_lts"
)

k_mapping <- c(
  "cd" = "k_cd",
  "lev" = "k_lev",
  "dfb" = "k_dfb",
  "mis_alpha" = "k_alpha",
  "mis_peel" = "k_peel_v2",
  "mis_sap" = "k_peel_sap",
  "mis_oracle" = "k_oracle"
)

overlap_mapping <- c(
  "cd" = "overlap_cd",
  "lev" = "overlap_lev",
  "dfb" = "overlap_dfb",
  "mis_alpha" = "overlap_mis_alpha",
  "mis_peel" = "overlap_peel_v2",
  "mis_sap" = "overlap_peel_sap",
  "mis_oracle" = "overlap_mis_oracle"
)

# ==============================================================================
# Main-paper exhibition sets
# ==============================================================================

# Estimation comparison shown in the main figures and tables.
MAIN_ESTIMATORS <- c(
  "full",
  "cd",
  "dfb",
  "mis_oracle",
  "mm",
  "lts"
)

# Detection methods shown in the main selected-k and overlap figures/tables.
MAIN_SELECTION_METHODS <- c(
  "cd",
  "dfb",
  "mis_oracle"
)

# ==============================================================================
# 4. Reshape iteration-level data
# ==============================================================================

coefficient_long <- reshape_mapped_columns(
  data = sim_primary,
  mapping = coef_mapping,
  value_name = "coefficient",
  metadata = estimator_meta,
  metadata_key = "estimator_id"
) %>%
  mutate(
    signed_error = coefficient - TRUE_BETA,
    squared_error = (coefficient - TRUE_BETA)^2
  )

bias_long <- reshape_mapped_columns(
  data = sim_primary,
  mapping = bias_mapping,
  value_name = "absolute_bias_recorded",
  metadata = estimator_meta,
  metadata_key = "estimator_id"
)

coverage_long <- reshape_mapped_columns(
  data = sim_primary,
  mapping = coverage_mapping,
  value_name = "coverage",
  metadata = estimator_meta,
  metadata_key = "estimator_id"
)

runtime_long <- reshape_mapped_columns(
  data = sim_primary,
  mapping = runtime_mapping,
  value_name = "runtime_seconds",
  metadata = estimator_meta,
  metadata_key = "estimator_id"
)

selection_long <- reshape_mapped_columns(
  data = sim_primary,
  mapping = k_mapping,
  value_name = "selected_k",
  metadata = selection_meta,
  metadata_key = "estimator_id"
)

selection_long <- selection_long %>%
  mutate(
    selected_prop = selected_k / n_obs
  )

overlap_long <- reshape_mapped_columns(
  data = sim_primary,
  mapping = overlap_mapping,
  value_name = "overlap",
  metadata = selection_meta,
  metadata_key = "estimator_id"
) %>%
  filter(outlier_method != "none")


# Validate that recorded absolute bias agrees with the coefficient columns.
bias_validation <- coefficient_long %>%
  select(
    iter,
    n_obs,
    design_k,
    contam_prop,
    x_type,
    error_type,
    outlier_method,
    estimator_id,
    coefficient
  ) %>%
  left_join(
    bias_long %>%
      select(
        iter,
        n_obs,
        design_k,
        contam_prop,
        x_type,
        error_type,
        outlier_method,
        estimator_id,
        absolute_bias_recorded
      ),
    by = c(
      "iter",
      "n_obs",
      "design_k",
      "contam_prop",
      "x_type",
      "error_type",
      "outlier_method",
      "estimator_id"
    )
  ) %>%
  mutate(
    absolute_bias_recomputed = abs(
      coefficient - TRUE_BETA
    ),
    difference =
      absolute_bias_recorded -
      absolute_bias_recomputed
  )

max_bias_difference <- if (
  any(is.finite(bias_validation$difference))
) {
  max(
    abs(bias_validation$difference),
    na.rm = TRUE
  )
} else {
  NA_real_
}


if (
  is.finite(max_bias_difference) &&
  max_bias_difference > 1e-8
) {
  warning(
    "Recorded bias columns differ from abs(coefficient - TRUE_BETA). ",
    "Maximum absolute difference: ",
    format(
      max_bias_difference,
      scientific = TRUE
    ),
    ". Publication summaries use coefficient-based recomputation."
  )
}


# ==============================================================================
# 5. Cell-level and equal-cell-weight broad summaries
# ==============================================================================

estimation_cell <- coefficient_long %>%
  group_by(
    n_obs,
    design_k,
    contam_prop,
    x_type,
    error_type,
    outlier_method,
    estimator_id,
    estimator_label,
    estimator_short,
    estimator_family,
    estimator_order
  ) %>%
  summarise(
    n_valid_coef = sum(is.finite(coefficient)),
    mean_coef = safe_mean(coefficient),
    sd_coef = safe_sd(coefficient),
    mcse_coef = safe_mcse(coefficient),
    
    mean_signed_bias = safe_mean(signed_error),
    mcse_signed_bias = safe_mcse(signed_error),
    
    mean_abs_bias = safe_mean(abs(signed_error)),
    sd_abs_bias = safe_sd(abs(signed_error)),
    mcse_abs_bias = safe_mcse(abs(signed_error)),
    
    mse = safe_mean(squared_error),
    mcse_mse = safe_mcse(squared_error),
    
    .groups = "drop"
  ) %>%
  mutate(
    rmse = sqrt(mse),
    mcse_rmse = ifelse(
      is.finite(rmse) & rmse > sqrt(.Machine$double.eps) &
        is.finite(mcse_mse),
      mcse_mse / (2 * rmse),
      NA_real_
    )
  )


coverage_cell <- coverage_long %>%
  group_by(
    n_obs,
    design_k,
    contam_prop,
    x_type,
    error_type,
    outlier_method,
    estimator_id
  ) %>%
  summarise(
    n_valid_coverage = sum(is.finite(coverage)),
    coverage_rate = safe_mean(coverage),
    mcse_coverage = safe_prop_mcse(coverage),
    .groups = "drop"
  )


estimation_cell <- estimation_cell %>%
  left_join(
    coverage_cell,
    by = c(
      "n_obs",
      "design_k",
      "contam_prop",
      "x_type",
      "error_type",
      "outlier_method",
      "estimator_id"
    )
  )


estimation_broad <- estimation_cell %>%
  group_by(
    outlier_method,
    estimator_id, estimator_label, estimator_short,
    estimator_family, estimator_order
  ) %>%
  summarise(
    n_design_cells = sum(is.finite(mean_coef)),
    
    mean_coef = safe_mean(mean_coef),
    mcse_coef = combined_cell_mcse(mcse_coef, mean_coef),
    
    mean_signed_bias = safe_mean(mean_signed_bias),
    mcse_signed_bias = combined_cell_mcse(
      mcse_signed_bias,
      mean_signed_bias
    ),
    
    mean_abs_bias = safe_mean(mean_abs_bias),
    mcse_abs_bias = combined_cell_mcse(
      mcse_abs_bias,
      mean_abs_bias
    ),
    
    mean_mse = safe_mean(mse),
    mcse_mse = combined_cell_mcse(mcse_mse, mse),
    
    coverage_rate = safe_mean(coverage_rate),
    mcse_coverage = combined_cell_mcse(
      mcse_coverage,
      coverage_rate
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    rmse = sqrt(mean_mse),
    mcse_rmse = ifelse(
      is.finite(rmse) & rmse > sqrt(.Machine$double.eps) &
        is.finite(mcse_mse),
      mcse_mse / (2 * rmse),
      NA_real_
    )
  )


selection_cell <- selection_long %>%
  group_by(
    n_obs,
    design_k,
    contam_prop,
    x_type,
    error_type,
    outlier_method,
    estimator_id,
    estimator_label,
    estimator_order,
    selection_order
  ) %>%
  summarise(
    mean_selected_prop = safe_mean(selected_prop),
    sd_selected_prop = safe_sd(selected_prop),
    mcse_selected_prop = safe_mcse(selected_prop),
    n_valid_k = sum(is.finite(selected_k)),
    mean_selected_k = safe_mean(selected_k),
    sd_selected_k = safe_sd(selected_k),
    mcse_selected_k = safe_mcse(selected_k),
    .groups = "drop"
  )


selection_broad <- selection_cell %>%
  group_by(
    outlier_method,
    estimator_id,
    estimator_label,
    estimator_order,
    selection_order
  ) %>%
  summarise(
    n_design_cells = sum(is.finite(mean_selected_k)),
    
    mean_selected_k = safe_mean(mean_selected_k),
    
    mcse_selected_k = combined_cell_mcse(
      mcse_selected_k,
      mean_selected_k
    ),
    
    mean_selected_prop = safe_mean(
      mean_selected_prop
    ),
    
    mcse_selected_prop = combined_cell_mcse(
      mcse_selected_prop,
      mean_selected_prop
    ),
    
    .groups = "drop"
  )

overlap_cell <- overlap_long %>%
  group_by(
    n_obs, design_k, contam_prop,
    x_type, error_type, outlier_method,
    estimator_id, estimator_label,
    estimator_order, selection_order
  ) %>%
  summarise(
    n_valid_overlap = sum(is.finite(overlap)),
    mean_overlap = safe_mean(overlap),
    sd_overlap = safe_sd(overlap),
    mcse_overlap = safe_mcse(overlap),
    .groups = "drop"
  )


overlap_broad <- overlap_cell %>%
  group_by(
    outlier_method,
    estimator_id, estimator_label,
    estimator_order, selection_order
  ) %>%
  summarise(
    n_design_cells = sum(is.finite(mean_overlap)),
    mean_overlap = safe_mean(mean_overlap),
    mcse_overlap = combined_cell_mcse(
      mcse_overlap,
      mean_overlap
    ),
    .groups = "drop"
  )


runtime_cell <- runtime_long %>%
  group_by(
    n_obs, design_k, contam_prop,
    x_type, error_type, outlier_method,
    estimator_id, estimator_label,
    estimator_order
  ) %>%
  summarise(
    n_valid_runtime = sum(is.finite(runtime_seconds)),
    mean_runtime = safe_mean(runtime_seconds),
    sd_runtime = safe_sd(runtime_seconds),
    mcse_runtime = safe_mcse(runtime_seconds),
    .groups = "drop"
  )


runtime_broad <- runtime_cell %>%
  group_by(
    outlier_method,
    estimator_id, estimator_label, estimator_order
  ) %>%
  summarise(
    n_design_cells = sum(is.finite(mean_runtime)),
    mean_runtime = safe_mean(mean_runtime),
    mcse_runtime = combined_cell_mcse(
      mcse_runtime,
      mean_runtime
    ),
    .groups = "drop"
  )


sap_cell <- sim_primary %>%
  group_by(
    n_obs,
    design_k,
    contam_prop,
    x_type,
    error_type,
    outlier_method
  ) %>%
  summarise(
    n_iter = n(),
    
    detection_rate = safe_mean(as.numeric(k_peel_sap > 0L)),
    mcse_detection = safe_prop_mcse(as.numeric(k_peel_sap > 0L)),
    
    mean_selected_prop = safe_mean( k_peel_sap / n_obs ),
    mcse_selected_prop = safe_mcse( k_peel_sap / n_obs ),
    
    mean_selected_k = safe_mean(k_peel_sap),
    mcse_selected_k = safe_mcse(k_peel_sap),
    
    exact_k_rate = safe_mean(as.numeric(k_peel_sap == set_size)),
    mcse_exact_k = safe_prop_mcse(as.numeric(k_peel_sap == set_size)),
    
    mean_overlap = safe_mean(overlap_peel_sap),
    mcse_overlap = safe_mcse(overlap_peel_sap),
    
    mean_iterations = safe_mean(peel_sap_iters),
    mcse_iterations = safe_mcse(peel_sap_iters),
    
    mean_final_p = safe_mean(peel_sap_final_p),
    mcse_final_p = safe_mcse(peel_sap_final_p),
    
    mean_min_p = safe_mean(peel_sap_min_p),
    mcse_min_p = safe_mcse(peel_sap_min_p),
    
    mean_peak_excess = safe_mean(peel_sap_peak_excess),
    mcse_peak_excess = safe_mcse(peel_sap_peak_excess),
    
    error_rate = safe_mean(as.numeric(peel_sap_stop == "error")),
    mcse_error = safe_prop_mcse(
      as.numeric(peel_sap_stop == "error")
    ),
    
    mean_abs_bias_sap = safe_mean(bias_mis_sap),
    mcse_abs_bias_sap = safe_mcse(bias_mis_sap),
    
    coverage_sap = safe_mean(cov_mis_sap),
    mcse_coverage_sap = safe_prop_mcse(cov_mis_sap),
    
    .groups = "drop"
  )


sap_broad <- sap_cell %>%
  group_by(outlier_method) %>%
  summarise(
    n_design_cells = n(),
    
    detection_rate = safe_mean(detection_rate),
    mcse_detection = combined_cell_mcse(
      mcse_detection,
      detection_rate
    ),
    
    mean_selected_prop = safe_mean(
      mean_selected_prop
    ),
    
    mcse_selected_prop = combined_cell_mcse(
      mcse_selected_prop,
      mean_selected_prop
    ),
    
    mean_selected_k = safe_mean(mean_selected_k),
    mcse_selected_k = combined_cell_mcse(
      mcse_selected_k,
      mean_selected_k
    ),
    
    exact_k_rate = safe_mean(exact_k_rate),
    mcse_exact_k = combined_cell_mcse(
      mcse_exact_k,
      exact_k_rate
    ),
    
    mean_overlap = safe_mean(mean_overlap),
    mcse_overlap = combined_cell_mcse(
      mcse_overlap,
      mean_overlap
    ),
    
    mean_iterations = safe_mean(mean_iterations),
    mcse_iterations = combined_cell_mcse(
      mcse_iterations,
      mean_iterations
    ),
    
    mean_final_p = safe_mean(mean_final_p),
    mcse_final_p = combined_cell_mcse(
      mcse_final_p,
      mean_final_p
    ),
    
    mean_min_p = safe_mean(mean_min_p),
    mcse_min_p = combined_cell_mcse(
      mcse_min_p,
      mean_min_p
    ),
    
    mean_peak_excess = safe_mean(mean_peak_excess),
    mcse_peak_excess = combined_cell_mcse(
      mcse_peak_excess,
      mean_peak_excess
    ),
    
    error_rate = safe_mean(error_rate),
    mcse_error = combined_cell_mcse(
      mcse_error,
      error_rate
    ),
    
    mean_abs_bias_sap = safe_mean(mean_abs_bias_sap),
    mcse_abs_bias_sap = combined_cell_mcse(
      mcse_abs_bias_sap,
      mean_abs_bias_sap
    ),
    
    coverage_sap = safe_mean(coverage_sap),
    mcse_coverage_sap = combined_cell_mcse(
      mcse_coverage_sap,
      coverage_sap
    ),
    
    .groups = "drop"
  )


# Save principal summaries before plotting.
utils::write.csv(
  estimation_cell,
  file.path(data_dir, "04_estimation_summary_by_cell.csv"),
  row.names = FALSE
)

utils::write.csv(
  estimation_broad,
  file.path(data_dir, "04_estimation_summary_equal_cell.csv"),
  row.names = FALSE
)

utils::write.csv(
  selection_cell,
  file.path(data_dir, "04_selected_k_summary_by_cell.csv"),
  row.names = FALSE
)

utils::write.csv(
  selection_broad,
  file.path(data_dir, "04_selected_k_summary_equal_cell.csv"),
  row.names = FALSE
)

utils::write.csv(
  overlap_cell,
  file.path(data_dir, "04_overlap_summary_by_cell.csv"),
  row.names = FALSE
)

utils::write.csv(
  overlap_broad,
  file.path(data_dir, "04_overlap_summary_equal_cell.csv"),
  row.names = FALSE
)

utils::write.csv(
  runtime_cell,
  file.path(data_dir, "04_runtime_summary_by_cell.csv"),
  row.names = FALSE
)

utils::write.csv(
  runtime_broad,
  file.path(data_dir, "04_runtime_summary_equal_cell.csv"),
  row.names = FALSE
)

utils::write.csv(
  sap_cell,
  file.path(data_dir, "04_sap_summary_by_cell.csv"),
  row.names = FALSE
)

utils::write.csv(
  sap_broad,
  file.path(data_dir, "04_sap_summary_equal_cell.csv"),
  row.names = FALSE
)


# ==============================================================================
# 6. Main Figure 1: coefficient-error distributions and heavy tails
# ==============================================================================

# Keep every finite coefficient error.
# No trimming, Winsorization, or quantile-based removal is applied.
coef_finite <- coefficient_long %>%
  filter(
    estimator_id %in% MAIN_ESTIMATORS,
    is.finite(signed_error)
  ) %>%
  select(
    iter,
    n_obs,
    design_k,
    contam_prop,
    x_type,
    error_type,
    outlier_method,
    estimator_id,
    estimator_label,
    estimator_order,
    signed_error
  ) %>%
  add_outlier_display() %>%
  mutate(
    estimator_label = factor(
      estimator_label,
      levels = estimator_meta$estimator_label
    )
  )


coef_distribution_summary <- coef_finite %>%
  group_by(
    outlier_method,
    outlier_label_plot,
    estimator_id,
    estimator_label,
    estimator_order
  ) %>%
  summarise(
    n_finite = n(),
    
    # Primary arithmetic-mean summary
    mean_error = safe_mean(signed_error),
    
    # Shown only to reveal skewness and mean-median separation
    median_error = stats::median(
      signed_error,
      na.rm = TRUE
    ),
    
    q05 = safe_quantile(signed_error, 0.05),
    q25 = safe_quantile(signed_error, 0.25),
    q75 = safe_quantile(signed_error, 0.75),
    q95 = safe_quantile(signed_error, 0.95),
    
    .groups = "drop"
  )


fig1_coef <- ggplot(
  coef_finite,
  aes(
    x = signed_error,
    y = estimator_label
  )
) +
  
  # Zero means that the estimated coefficient equals TRUE_BETA
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    colour = COL_ORANGE
  ) +
  
  # Central 90% interval
  geom_segment(
    data = coef_distribution_summary,
    aes(
      x = q05,
      xend = q95,
      y = estimator_label,
      yend = estimator_label
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    colour = COL_GREY
  ) +
  
  # Interquartile interval
  geom_segment(
    data = coef_distribution_summary,
    aes(
      x = q25,
      xend = q75,
      y = estimator_label,
      yend = estimator_label
    ),
    inherit.aes = FALSE,
    linewidth = 2.2,
    colour = COL_BLUE_LIGHT
  ) +
  
  # Iteration-level coefficient errors
  geom_point(
    aes(shape = "Monte Carlo draw"),
    position = position_jitter(
      width = 0,
      height = 0.11,
      seed = 84
    ),
    alpha = 0.12,
    size = 0.40,
    colour = COL_BLUE
  ) +
  
  # Arithmetic mean across Monte Carlo draws
  geom_point(
    data = coef_distribution_summary,
    aes(
      x = mean_error,
      y = estimator_label,
      shape = "Mean"
    ),
    inherit.aes = FALSE,
    size = 2.4,
    stroke = 0.5,
    fill = COL_BLUE_DARK,
    colour = COL_BLACK
  ) +
  
  # Median across Monte Carlo draws
  geom_point(
    data = coef_distribution_summary,
    aes(
      x = median_error,
      y = estimator_label,
      shape = "Median"
    ),
    inherit.aes = FALSE,
    size = 2.3,
    stroke = 0.6,
    fill = "white",
    colour = COL_BLACK
  ) +
  
  facet_wrap(
    ~ outlier_label_plot,
    ncol = 2,
    scales = "fixed"
  ) +
  
  # Preserves negative values and zero while compressing extreme tails
  scale_x_continuous(
    trans = scales::pseudo_log_trans(
      base = 10,
      sigma = 0.05
    ),
    breaks = c(
      -3, -1, -0.3, -0.1, -0.03,
      0,
      0.03, 0.1, 0.3, 1, 3
    ),
    labels = c(
      "-3", "-1", "-0.3", "-0.1", "-0.03",
      "0",
      "0.03", "0.1", "0.3", "1", "3"
    )
  ) +
  
  scale_shape_manual(
    name = NULL,
    breaks = c(
      "Monte Carlo draw",
      "Mean",
      "Median"
    ),
    values = c(
      "Monte Carlo draw" = 16,
      "Mean" = 21,
      "Median" = 23
    )
  ) +
  
  guides(
    shape = guide_legend(
      direction = "horizontal",
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        alpha = c(0.30, 1, 1),
        size = c(1.7, 3.0, 3.0),
        colour = c(
          COL_BLUE,
          COL_BLACK,
          COL_BLACK
        ),
        fill = c(
          NA,
          COL_BLUE_DARK,
          "white"
        ),
        stroke = c(0, 0.5, 0.6)
      )
    )
  ) +
  
  labs(
    x = expression(
      hat(beta) - beta[0] ~ "(signed pseudo-log scale)"
    ),
    y = NULL
  ) +
  
  theme_distribution(base_size = 9.5) +
  
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 8.5),
    legend.margin = margin(b = 3),
    legend.box.margin = margin(b = 2),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_text(
      size = 8,
      margin = margin(t = 4)
    ),
    
    axis.title.x = element_text(
      margin = margin(t = 10)
    ),
    
    plot.margin = margin(
      t = 8,
      r = 16,
      b = 10,
      l = 8
    )
  )


save_plot(
  fig1_coef,
  file.path(
    fig_main_dir,
    "04_fig1_coefficient_error_distributions.pdf"
  ),
  width = 12.0,
  height = 8.0
)


utils::write.csv(
  coef_distribution_summary,
  file.path(
    data_dir,
    "04_fig1_coe_distribution_summary.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 7. Main Figure 2: DGP conditions associated with large MIS-oracle errors
# ==============================================================================

# ------------------------------------------------------------------------------
# Figure 2 statistical settings
# ------------------------------------------------------------------------------

# Since TRUE_BETA = 1, an absolute coefficient error greater than 0.5
# represents an error of at least half the true coefficient magnitude.
TAIL_ERROR_THRESHOLD <- 0.5

# Used to identify problematic DGP cells for the supplementary
# all-estimator heatmap.
PROBLEM_RATE_CUTOFF <- 0.05


# ------------------------------------------------------------------------------
# Figure 2 adjustable graphical settings
# ------------------------------------------------------------------------------

FIG2_MAIN_NCOL <- 2L

FIG2_MAIN_TILE_WIDTH  <- 0.96
FIG2_MAIN_TILE_HEIGHT <- 0.92
FIG2_MAIN_TILE_TEXT_SIZE <- 3.0

FIG2_MAIN_X_TEXT_SIZE <- 8.2
FIG2_MAIN_Y_TEXT_SIZE <- 8.8
FIG2_MAIN_STRIP_TEXT_SIZE <- 10.5

FIG2_MAIN_LEGEND_WIDTH_CM  <- 9.5
FIG2_MAIN_LEGEND_HEIGHT_CM <- 0.45

FIG2_MAIN_WIDTH_IN  <- 11.2
FIG2_MAIN_HEIGHT_IN <- 7.4


# ------------------------------------------------------------------------------
# Remove obsolete Figure 2 outputs
# ------------------------------------------------------------------------------

obsolete_fig2_files <- c(
  file.path(
    fig_main_dir,
    "04_fig2_bias_coverage_tradeoff.pdf"
  ),
  file.path(
    data_dir,
    "04_fig2_bias_coverage_tradeoff_data.csv"
  ),
  file.path(
    fig_main_dir,
    "04_fig2_large_error_dgp_heatmap.pdf"
  ),
  file.path(
    data_dir,
    "04_fig2_large_error_dgp_heatmap_data.csv"
  )
)

obsolete_fig2_files <- obsolete_fig2_files[
  file.exists(obsolete_fig2_files)
]

if (length(obsolete_fig2_files) > 0L) {
  unlink(obsolete_fig2_files)
  
  cat(
    "Removed obsolete Figure 2 file(s):\n",
    paste0(
      "  - ",
      obsolete_fig2_files,
      collapse = "\n"
    ),
    "\n",
    sep = ""
  )
}


# ------------------------------------------------------------------------------
# Summarise large-error frequency for every estimator and exact DGP cell
# ------------------------------------------------------------------------------

coefficient_long_all <- reshape_mapped_columns(
  data = sim_all,
  mapping = coef_mapping,
  value_name = "coefficient",
  metadata = estimator_meta,
  metadata_key = "estimator_id"
) %>%
  mutate(
    signed_error = coefficient - TRUE_BETA,
    squared_error = (coefficient - TRUE_BETA)^2
  )

tail_cell_summary <- coefficient_long_all %>%
  mutate(
    finite_error = is.finite(signed_error),
    
    large_error = ifelse(
      finite_error,
      abs(signed_error) > TAIL_ERROR_THRESHOLD,
      NA
    ),
    
    positive_large_error = ifelse(
      finite_error,
      signed_error > TAIL_ERROR_THRESHOLD,
      NA
    ),
    
    negative_large_error = ifelse(
      finite_error,
      signed_error < -TAIL_ERROR_THRESHOLD,
      NA
    )
  ) %>%
  group_by(
    n_obs,
    design_k,
    contam_prop,
    x_type,
    error_type,
    outlier_method,
    estimator_id,
    estimator_label,
    estimator_order
  ) %>%
  summarise(
    n_total = n(),
    n_finite = sum(finite_error),
    
    exceedance_rate = safe_mean(
      as.numeric(large_error)
    ),
    
    mcse_exceedance = safe_prop_mcse(
      as.numeric(large_error)
    ),
    
    positive_exceedance_rate = safe_mean(
      as.numeric(positive_large_error)
    ),
    
    negative_exceedance_rate = safe_mean(
      as.numeric(negative_large_error)
    ),
    
    nonfinite_rate = mean(!finite_error),
    
    .groups = "drop"
  )

tail_equal_grid_summary <- tail_cell_summary %>%
  group_by(
    x_type,
    error_type,
    outlier_method,
    estimator_id,
    estimator_label,
    estimator_order
  ) %>%
  summarise(
    exceedance_rate = safe_mean(
      exceedance_rate
    ),
    
    mcse_exceedance = combined_cell_mcse(
      mcse_exceedance,
      exceedance_rate
    ),
    
    positive_exceedance_rate = safe_mean(
      positive_exceedance_rate
    ),
    
    negative_exceedance_rate = safe_mean(
      negative_exceedance_rate
    ),
    
    nonfinite_rate = safe_mean(
      nonfinite_rate
    ),
    
    .groups = "drop"
  )


# ------------------------------------------------------------------------------
# Identify problematic DGP cells for the supplementary all-estimator figure
# ------------------------------------------------------------------------------

problem_dgp_cells <- tail_equal_grid_summary %>%
  group_by(
    x_type,
    error_type,
    outlier_method
  ) %>%
  summarise(
    maximum_exceedance_rate = if (
      any(is.finite(exceedance_rate))
    ) {
      max(
        exceedance_rate,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    
    .groups = "drop"
  ) %>%
  filter(
    is.finite(maximum_exceedance_rate),
    maximum_exceedance_rate >= PROBLEM_RATE_CUTOFF
  )


# ------------------------------------------------------------------------------
# Shared cell-label function
# ------------------------------------------------------------------------------

add_tail_annotations <- function(data) {
  data %>%
    mutate(
      dominant_tail = case_when(
        !is.finite(exceedance_rate) ~ "",
        
        exceedance_rate <= 0 ~ "",
        
        positive_exceedance_rate -
          negative_exceedance_rate >= 0.01 ~ "+",
        
        negative_exceedance_rate -
          positive_exceedance_rate >= 0.01 ~ "-",
        
        TRUE ~ "\u00B1"
      ),
      
      tile_label = case_when(
        !is.finite(exceedance_rate) ~ "NA",
        
        dominant_tail == "" ~ scales::percent(
          exceedance_rate,
          accuracy = 1
        ),
        
        TRUE ~ paste0(
          scales::percent(
            exceedance_rate,
            accuracy = 1
          ),
          "\n",
          dominant_tail
        )
      ),
      
      label_colour = ifelse(
        is.finite(exceedance_rate) &
          exceedance_rate >= 0.25,
        "white",
        COL_BLACK
      )
    )
}


# ------------------------------------------------------------------------------
# Shared blue-white-orange fill scale
# ------------------------------------------------------------------------------

tail_fill_scale <- function(legend_title) {
  scale_fill_gradientn(
    colours = c(
      COL_BLUE_LIGHT,
      COL_NEUTRAL,
      COL_ORANGE,
      COL_ORANGE_DARK
    ),
    
    values = scales::rescale(
      c(
        0,
        PROBLEM_RATE_CUTOFF,
        0.25,
        1
      )
    ),
    
    limits = c(0, 1),
    oob = scales::squish,
    na.value = COL_GREY_LIGHT,
    
    # Avoid placing 0%, 5%, and 10% too close together.
    breaks = c(
      0,
      0.10,
      0.25,
      0.50,
      1
    ),
    
    labels = scales::label_percent(
      accuracy = 1
    ),
    
    name = legend_title
  )
}


# ------------------------------------------------------------------------------
# MIS-oracle data for the main-paper heatmap
# ------------------------------------------------------------------------------

tail_mis_oracle_data <- tail_equal_grid_summary %>%
  filter(
    estimator_id == "mis_oracle"
  ) %>%
  add_outlier_display() %>%
  mutate(
    predictor_label = factor(
      unname(
        x_labels_table[
          as.character(x_type)
        ]
      ),
      levels = rev(
        unname(
          x_labels_table[x_order]
        )
      )
    ),
    
    error_label = factor(
      unname(
        error_labels_table[
          as.character(error_type)
        ]
      ),
      levels = unname(
        error_labels_table[error_order]
      )
    )
  ) %>%
  add_tail_annotations() %>%
  arrange(
    outlier_method,
    x_type,
    error_type
  )


# ------------------------------------------------------------------------------
# Main Figure 2: MIS-oracle only
# ------------------------------------------------------------------------------

fig2_mis_oracle_tail <- ggplot(
  tail_mis_oracle_data,
  aes(
    x = error_label,
    y = predictor_label,
    fill = exceedance_rate
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.55,
    width = FIG2_MAIN_TILE_WIDTH,
    height = FIG2_MAIN_TILE_HEIGHT
  ) +
  
  geom_text(
    aes(
      label = tile_label,
      colour = label_colour
    ),
    size = FIG2_MAIN_TILE_TEXT_SIZE,
    lineheight = 0.88,
    show.legend = FALSE
  ) +
  
  facet_wrap(
    ~outlier_label_plot,
    ncol = FIG2_MAIN_NCOL,
    drop = TRUE
  ) +
  
  tail_fill_scale(
    paste0(
      "MIS-oracle probability that absolute\n",
      "coefficient error exceeds ",
      TAIL_ERROR_THRESHOLD
    )
  ) +
  
  scale_colour_identity() +
  
  scale_x_discrete(
    drop = FALSE,
    expand = expansion(
      mult = c(0.01, 0.01)
    )
  ) +
  
  scale_y_discrete(
    drop = FALSE,
    expand = expansion(
      mult = c(0.02, 0.02)
    )
  ) +
  
  labs(
    x = "Error distribution",
    y = "Predictor distribution"
  ) +
  
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      barwidth = grid::unit(
        FIG2_MAIN_LEGEND_WIDTH_CM,
        "cm"
      ),
      barheight = grid::unit(
        FIG2_MAIN_LEGEND_HEIGHT_CM,
        "cm"
      ),
      ticks = TRUE
    )
  ) +
  
  theme_heatmap(base_size = 10) +
  
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1,
      size = FIG2_MAIN_X_TEXT_SIZE,
      colour = COL_BLACK
    ),
    
    axis.text.y = element_text(
      size = FIG2_MAIN_Y_TEXT_SIZE,
      colour = COL_BLACK
    ),
    
    axis.title.x = element_text(
      colour = COL_BLACK,
      margin = margin(t = 10)
    ),
    
    axis.title.y = element_text(
      colour = COL_BLACK,
      margin = margin(r = 10)
    ),
    
    strip.text = element_text(
      size = FIG2_MAIN_STRIP_TEXT_SIZE,
      face = "bold",
      colour = COL_BLACK,
      margin = margin(
        t = 5,
        b = 5
      )
    ),
    
    legend.title = element_text(
      size = 9.5,
      colour = COL_BLACK,
      hjust = 0.5
    ),
    
    legend.text = element_text(
      size = 8.5,
      colour = COL_BLACK
    ),
    
    panel.spacing = grid::unit(
      1.2,
      "lines"
    ),
    
    plot.margin = margin(
      t = 8,
      r = 14,
      b = 14,
      l = 10
    )
  )


save_plot(
  fig2_mis_oracle_tail,
  file.path(
    fig_main_dir,
    "04_fig2_mis_oracle_large_error_dgp_heatmap.pdf"
  ),
  width = FIG2_MAIN_WIDTH_IN,
  height = FIG2_MAIN_HEIGHT_IN
)


# ------------------------------------------------------------------------------
# Save Figure 2 data
# ------------------------------------------------------------------------------

utils::write.csv(
  tail_mis_oracle_data,
  file.path(
    data_dir,
    "04_fig2_mis_oracle_large_error_dgp_heatmap_data.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  tail_cell_summary,
  file.path(
    data_dir,
    "04_fig2_large_error_summary_all_cells.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  problem_dgp_cells,
  file.path(
    data_dir,
    "04_fig2_problem_dgp_cells.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 8. Main Figure 3: detection recall and precision distributions
# ==============================================================================

detection_join_keys <- c(
  "iter",
  "n_obs",
  "design_k",
  "contam_prop",
  "x_type",
  "error_type",
  "outlier_method",
  "estimator_id"
)


# Join selected-k and overlap information at the Monte Carlo iteration level.
detection_quality_iteration <- selection_long %>%
  filter(
    outlier_method != "none",
    estimator_id %in% MAIN_SELECTION_METHODS,
    is.finite(selected_k)
  ) %>%
  select(
    all_of(detection_join_keys),
    estimator_label,
    estimator_order,
    selection_order,
    selected_k
  ) %>%
  inner_join(
    overlap_long %>%
      filter(
        estimator_id %in% MAIN_SELECTION_METHODS,
        is.finite(overlap)
      ) %>%
      select(
        all_of(detection_join_keys),
        overlap
      ),
    by = detection_join_keys
  ) %>%
  mutate(
    # overlap is the fraction of the true injected coalition recovered.
    recall = pmin(
      1,
      pmax(0, overlap)
    ),
    
    true_positive =
      recall * design_k,
    
    # Precision measures how much of the selected deletion set
    # actually belongs to the injected coalition.
    precision = ifelse(
      selected_k > 0,
      pmin(
        1,
        pmax(
          0,
          true_positive / selected_k
        )
      ),
      NA_real_
    )
  )


# DGP-cell summaries.
detection_quality_cell <- detection_quality_iteration %>%
  group_by(
    n_obs,
    design_k,
    contam_prop,
    x_type,
    error_type,
    outlier_method,
    estimator_id,
    estimator_label,
    estimator_order,
    selection_order
  ) %>%
  summarise(
    n_valid_selected_k =
      sum(is.finite(selected_k)),
    
    n_valid_recall =
      sum(is.finite(recall)),
    
    n_valid_precision =
      sum(is.finite(precision)),
    
    mean_selected_k =
      safe_mean(selected_k),
    
    sd_selected_k =
      safe_sd(selected_k),
    
    mcse_selected_k =
      safe_mcse(selected_k),
    
    mean_recall =
      safe_mean(recall),
    
    sd_recall =
      safe_sd(recall),
    
    mcse_recall =
      safe_mcse(recall),
    
    mean_precision =
      safe_mean(precision),
    
    sd_precision =
      safe_sd(precision),
    
    mcse_precision =
      safe_mcse(precision),
    
    .groups = "drop"
  )


# Equal-cell-weight broad summaries used by Main Table 3.
detection_quality_broad <- detection_quality_cell %>%
  group_by(
    outlier_method,
    estimator_id,
    estimator_label,
    estimator_order,
    selection_order
  ) %>%
  summarise(
    n_design_cells =
      sum(is.finite(mean_recall)),
    
    broad_mean_selected_k =
      safe_mean(mean_selected_k),
    
    broad_mcse_selected_k =
      combined_cell_mcse(
        mcse_selected_k,
        mean_selected_k
      ),
    
    broad_mean_recall =
      safe_mean(mean_recall),
    
    broad_mcse_recall =
      combined_cell_mcse(
        mcse_recall,
        mean_recall
      ),
    
    broad_mean_precision =
      safe_mean(mean_precision),
    
    broad_mcse_precision =
      combined_cell_mcse(
        mcse_precision,
        mean_precision
      ),
    
    .groups = "drop"
  ) %>%
  transmute(
    outlier_method,
    estimator_id,
    estimator_label,
    estimator_order,
    selection_order,
    n_design_cells,
    
    mean_selected_k =
      broad_mean_selected_k,
    
    mcse_selected_k =
      broad_mcse_selected_k,
    
    mean_recall =
      broad_mean_recall,
    
    mcse_recall =
      broad_mcse_recall,
    
    mean_precision =
      broad_mean_precision,
    
    mcse_precision =
      broad_mcse_precision
  )


# Convert recall and precision into one plotting variable.
detection_quality_plot_data <- detection_quality_cell %>%
  select(
    n_obs,
    design_k,
    contam_prop,
    x_type,
    error_type,
    outlier_method,
    estimator_id,
    estimator_label,
    estimator_order,
    selection_order,
    mean_recall,
    mean_precision
  ) %>%
  pivot_longer(
    cols = c(
      mean_recall,
      mean_precision
    ),
    names_to = "metric_id",
    values_to = "metric_value"
  ) %>%
  filter(
    is.finite(metric_value)
  ) %>%
  mutate(
    metric_label = factor(
      recode(
        metric_id,
        "mean_recall" = "Recall",
        "mean_precision" = "Precision"
      ),
      levels = c(
        "Recall",
        "Precision"
      )
    )
  ) %>%
  add_outlier_display() %>%
  mutate(
    estimator_label = factor(
      estimator_label,
      levels = selection_meta$estimator_label[
        match(
          MAIN_SELECTION_METHODS,
          selection_meta$estimator_id
        )
      ]
    )
  )


detection_quality_distribution_summary <-
  detection_quality_plot_data %>%
  group_by(
    outlier_method,
    outlier_label_plot,
    estimator_id,
    estimator_label,
    estimator_order,
    selection_order,
    metric_id,
    metric_label
  ) %>%
  summarise(
    n_valid =
      sum(is.finite(metric_value)),
    
    mean_metric =
      safe_mean(metric_value),
    
    sd_metric =
      safe_sd(metric_value),
    
    .groups = "drop"
  )


fig3_detection_quality <- ggplot(
  detection_quality_plot_data,
  aes(
    x = metric_value,
    y = estimator_label
  )
) +
  
  geom_violin(
    fill = COL_VIOLIN,
    colour = COL_GREY_DARK,
    linewidth = 0.25,
    scale = "width",
    trim = TRUE,
    na.rm = TRUE
  ) +
  
  geom_segment(
    data = detection_quality_distribution_summary,
    aes(
      x = pmax(
        0,
        mean_metric - sd_metric
      ),
      xend = pmin(
        1,
        mean_metric + sd_metric
      ),
      y = estimator_label,
      yend = estimator_label
    ),
    inherit.aes = FALSE,
    colour = COL_BLUE,
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  
  geom_point(
    data = detection_quality_distribution_summary,
    aes(
      x = mean_metric,
      y = estimator_label
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.2,
    stroke = 0.45,
    fill = COL_BLUE,
    colour = COL_BLACK,
    na.rm = TRUE
  ) +
  
  geom_vline(
    xintercept = 1,
    colour = COL_ORANGE,
    linetype = "dashed",
    linewidth = 0.75
  ) +
  
  facet_grid(
    metric_label ~ outlier_label_plot
  ) +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(
      0,
      1,
      by = 0.2
    ),
    labels =
      scales::label_percent(
        accuracy = 1
      )
  ) +
  
  labs(
    x = "Proportion",
    y = NULL
  ) +
  
  theme_distribution() +
  
  theme(
    legend.position = "none"
  )


save_plot(
  fig3_detection_quality,
  file.path(
    fig_main_dir,
    "04_fig3_detection_quality_distributions.pdf"
  ),
  width = 11.2,
  height = 7.2
)


utils::write.csv(
  detection_quality_distribution_summary,
  file.path(
    data_dir,
    "04_fig3_detection_quality_distribution_summary.csv"
  ),
  row.names = FALSE
)


utils::write.csv(
  detection_quality_cell,
  file.path(
    data_dir,
    "04_detection_quality_summary_by_cell.csv"
  ),
  row.names = FALSE
)


utils::write.csv(
  detection_quality_broad,
  file.path(
    data_dir,
    "04_detection_quality_summary_equal_cell.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 10. Main tables
# ==============================================================================

estimation_broad_main <- estimation_broad %>%
  filter(
    estimator_id %in% MAIN_ESTIMATORS
  )

# Table 1: arithmetic mean coefficients.
tab1_mean_coef <- make_metric_wide_table(
  data = estimation_broad_main,
  mean_column = "mean_coef",
  mcse_column = "mcse_coef",
  formatter = function(x, se) fmt_mean_mcse(x, se, digits = 3L),
  include_clean = TRUE
)

write_tex_table(
  data = tab1_mean_coef,
  tex_path = file.path(
    tab_main_dir,
    "04_tab1_mean_coefficients.tex"
  ),
  caption = paste0(
    "Arithmetic mean coefficient estimates by estimator and contamination ",
    "mechanism. The true coefficient is beta0 = 1. Results give equal ",
    "weight to each sample-size, contamination-proportion, ",
    "predictor-distribution, and error-distribution design cell. ",
    "Monte Carlo standard errors are in parentheses."
  ),
  label = "tab:robust-mean-coefficients",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrr"
)


# Table 2a: mean absolute bias.
tab2a_bias <- make_metric_wide_table(
  data = estimation_broad_main,
  mean_column = "mean_abs_bias",
  mcse_column = "mcse_abs_bias",
  formatter = function(x, se) fmt_mean_mcse(x, se, digits = 3L),
  include_clean = TRUE
)

write_tex_table(
  data = tab2a_bias,
  tex_path = file.path(
    tab_main_dir,
    "04_tab2a_mean_absolute_bias.tex"
  ),
  caption = paste0(
    "Mean absolute coefficient bias by estimator and contamination mechanism. ",
    "Results give equal weight to each sample-size, ",
    "contamination-proportion, predictor-distribution, and ",
    "error-distribution design cell. Monte Carlo standard errors are ",
    "in parentheses."
  ),
  label = "tab:robust-mean-absolute-bias",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrr"
)


# Table 2b: RMSE.
tab2b_rmse <- make_metric_wide_table(
  data = estimation_broad_main,
  mean_column = "rmse",
  mcse_column = "mcse_rmse",
  formatter = function(x, se) fmt_mean_mcse(x, se, digits = 3L),
  include_clean = TRUE
)

write_tex_table(
  data = tab2b_rmse,
  tex_path = file.path(
    tab_main_dir,
    "04_tab2b_rmse.tex"
  ),
  caption = paste0(
    "Root mean squared error by estimator and contamination mechanism. ",
    "Results give equal weight to each sample-size, ",
    "contamination-proportion, predictor-distribution, and ",
    "error-distribution design cell. Delta-method Monte Carlo standard ",
    "errors are in parentheses."
  ),
  label = "tab:robust-rmse",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrr"
)


# Table 2c: empirical coverage.
tab2c_coverage <- make_metric_wide_table(
  data = estimation_broad_main,
  mean_column = "coverage_rate",
  mcse_column = "mcse_coverage",
  formatter = function(x, se) fmt_pct_mcse(x, se, digits = 1L),
  include_clean = TRUE
)

write_tex_table(
  data = tab2c_coverage,
  tex_path = file.path(
    tab_main_dir,
    "04_tab2c_coverage.tex"
  ),
  caption = paste0(
    "Empirical coverage of nominal 95 percent confidence intervals. ",
    "Results give equal weight to each sample-size, ",
    "contamination-proportion, predictor-distribution, and ",
    "error-distribution design cell. Monte Carlo standard errors in ",
    "percentage points are in parentheses."
  ),
  label = "tab:robust-coverage",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrr",
  escape_cells = FALSE
)

# Table 3: detection efficiency.
tab3_detection_efficiency <- detection_quality_broad %>%
  add_outlier_display() %>%
  arrange(
    outlier_label_table,
    selection_order
  ) %>%
  transmute(
    Scenario =
      as.character(outlier_label_table),
    
    Method =
      as.character(estimator_label),
    
    `Mean selected k` =
      fmt_mean_mcse(
        mean_selected_k,
        mcse_selected_k,
        digits = 1L
      ),
    
    Recall =
      fmt_pct_mcse(
        mean_recall,
        mcse_recall,
        digits = 1L
      ),
    
    Precision =
      fmt_pct_mcse(
        mean_precision,
        mcse_precision,
        digits = 1L
      )
  )


write_tex_table(
  data = tab3_detection_efficiency,
  tex_path = file.path(
    tab_main_dir,
    "04_tab3_detection_efficiency.tex"
  ),
  caption = paste0(
    "Detection-set efficiency by contamination mechanism. ",
    "Mean selected k is the mean number of observations removed; ",
    "recall is the fraction of the injected coalition recovered, and ",
    "precision is the fraction of selected observations belonging to the ",
    "injected coalition. The contaminated designs use N = 5000 and true ",
    "coalition size k = 50 (1 percent contamination). MIS with oracle k is ",
    "supplied the true coalition size, so its selected-k entry represents ",
    "the fixed deletion budget rather than an estimated coalition size. ",
    "Results give equal weight to predictor-distribution and ",
    "error-distribution design cells. Monte Carlo standard errors are ",
    "in parentheses."
  ),
  label = "tab:robust-detection-efficiency",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "llrrr",
  escape_cells = FALSE
)

# ==============================================================================
# 11. Supplementary Figure A1: absolute-bias distributions
# ==============================================================================

bias_plot_data <- coefficient_long %>%
  mutate(absolute_bias = abs(coefficient - TRUE_BETA)) %>%
  filter(
    estimator_id %in% MAIN_ESTIMATORS,
    is.finite(absolute_bias)
  ) %>%
  add_outlier_display() %>%
  mutate(
    estimator_label = factor(
      estimator_label,
      levels = estimator_meta$estimator_label
    )
  )


bias_distribution_summary <- bias_plot_data %>%
  group_by(
    outlier_method,
    outlier_label_plot,
    estimator_id, estimator_label, estimator_order
  ) %>%
  summarise(
    n_valid = n(),
    mean_absolute_bias = safe_mean(absolute_bias),
    sd_absolute_bias = safe_sd(absolute_bias),
    .groups = "drop"
  )


figA1_bias <- ggplot(
  bias_plot_data,
  aes(x = absolute_bias, y = estimator_label)
) +
  geom_violin(
    fill = COL_VIOLIN,
    colour = COL_GREY_DARK,
    linewidth = 0.25,
    scale = "width",
    trim = TRUE
  ) +
  geom_segment(
    data = bias_distribution_summary,
    aes(
      x = pmax(0, mean_absolute_bias - sd_absolute_bias),
      xend = mean_absolute_bias + sd_absolute_bias,
      y = estimator_label,
      yend = estimator_label
    ),
    inherit.aes = FALSE,
    colour = COL_BLUE,
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  geom_point(
    data = bias_distribution_summary,
    aes(
      x = mean_absolute_bias,
      y = estimator_label
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.2,
    stroke = 0.45,
    fill = COL_BLUE,
    colour = COL_BLACK,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ outlier_label_plot,
    ncol = 2
  ) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(
      base = 10,
      sigma = 0.001
    ),
    
    # Short, interpretable labels for the very wide bias range
    breaks = c(
      0,
      1,
      1e2,
      1e4,
      1e6,
      1e8
    ),
    
    labels = c(
      "0",
      "1",
      "100",
      "10K",
      "1M",
      "100M"
    ),
    
    expand = expansion(
      mult = c(0.01, 0.03)
    ),
    
    # Final protection against accidental overlap
    guide = guide_axis(
      check.overlap = TRUE
    )
  ) +
  labs(
    x = "Absolute coefficient bias",
    y = NULL
  ) +
  theme_distribution() +
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      size = 8,
      colour = COL_BLACK,
      margin = margin(t = 5)
    ),
    
    axis.title.x = element_text(
      colour = COL_BLACK,
      margin = margin(t = 10)
    ),
    
    axis.text.y = element_text(
      colour = COL_BLACK
    ),
    
    # More separation between the two facet columns
    panel.spacing.x = grid::unit(
      1.6,
      "lines"
    ),
    
    plot.margin = margin(
      t = 8,
      r = 18,
      b = 16,
      l = 8
    )
  )


save_plot(
  figA1_bias,
  file.path(
    fig_supp_dir,
    "04_figA1_absolute_bias_distributions.pdf"
  ),
  width = 11.5,
  height = 7.6
)

utils::write.csv(
  bias_distribution_summary,
  file.path(data_dir, "04_figA1_absolute_bias_summary.csv"),
  row.names = FALSE
)


# ==============================================================================
# 12. Supplementary Figure A2: mean bias by error distribution
# ==============================================================================

selected_methods_for_sensitivity <- MAIN_ESTIMATORS

bias_by_error <- estimation_cell %>%
  filter(estimator_id %in% selected_methods_for_sensitivity) %>%
  group_by(
    error_type, outlier_method,
    estimator_id, estimator_label, estimator_order
  ) %>%
  summarise(
    mean_abs_bias = safe_mean(mean_abs_bias),
    .groups = "drop"
  ) %>%
  add_outlier_display() %>%
  mutate(
    error_label = factor(
      unname(error_labels_plot[as.character(error_type)]),
      levels = unname(error_labels_plot[error_order])
    ),
    estimator_id = factor(
      estimator_id,
      levels = estimator_meta$estimator_id
    )
  )

compact_pseudolog_breaks <- function(limits) {
  finite_limits <- limits[is.finite(limits)]
  
  if (length(finite_limits) == 0L) {
    return(0)
  }
  
  upper <- max(finite_limits)
  
  if (!is.finite(upper) || upper <= 0) {
    return(0)
  }
  
  exponent_max <- ceiling(log10(upper))
  exponent_min <- floor(
    log10(max(upper / 1000, 1e-12))
  )
  
  exponents <- seq(
    exponent_min,
    exponent_max
  )
  
  if (length(exponents) > 4L) {
    exponents <- exponents[
      unique(round(seq(
        1,
        length(exponents),
        length.out = 4
      )))
    ]
  }
  
  sort(unique(c(
    0,
    10^exponents
  )))
}


compact_scientific_labels <- function(x) {
  ifelse(
    x == 0,
    "0",
    formatC(
      x,
      format = "e",
      digits = 0
    )
  )
}

figA2_bias_error <- ggplot(
  bias_by_error,
  aes(
    x = error_label,
    y = mean_abs_bias,
    group = estimator_id,
    colour = estimator_id,
    shape = estimator_id,
    linetype = estimator_id
  )
) +
  geom_line(
    linewidth = 0.55,
    na.rm = TRUE
  ) +
  geom_point(
    size = 2.0,
    stroke = 0.6,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ outlier_label_plot,
    ncol = 2,
    scales = "free_y"
  ) +
  scale_x_discrete(
    drop = FALSE,
    guide = guide_axis(
      n.dodge = 2,
      check.overlap = FALSE
    )
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(
      base = 10,
      sigma = 0.001
    ),
    breaks = compact_pseudolog_breaks,
    labels = compact_scientific_labels,
    expand = expansion(
      mult = c(0.02, 0.08)
    )
  ) +
  scale_colour_manual(
    values = method_colors,
    breaks = MAIN_ESTIMATORS,
    labels = estimator_meta$estimator_label[
      match(
        MAIN_ESTIMATORS,
        estimator_meta$estimator_id
      )
    ],
    drop = TRUE
  ) +
  scale_shape_manual(
    values = method_shapes,
    breaks = MAIN_ESTIMATORS,
    labels = estimator_meta$estimator_label[
      match(
        MAIN_ESTIMATORS,
        estimator_meta$estimator_id
      )
    ],
    drop = TRUE
  ) +
  scale_linetype_manual(
    values = method_linetypes,
    breaks = MAIN_ESTIMATORS,
    labels = estimator_meta$estimator_label[
      match(
        MAIN_ESTIMATORS,
        estimator_meta$estimator_id
      )
    ],
    drop = TRUE
  ) +
  labs(
    x = "Error distribution",
    y = "Mean absolute bias",
    colour = "Estimator",
    shape = "Estimator",
    linetype = "Estimator"
  ) +
  theme_paper() +
  theme_paper() +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 1,
      size = 8,
      lineheight = 0.9,
      colour = COL_BLACK,
      margin = margin(t = 5)
    ),
    
    axis.text.y = element_text(
      size = 8,
      colour = COL_BLACK
    ),
    
    axis.title.x = element_text(
      colour = COL_BLACK,
      margin = margin(t = 12)
    ),
    
    axis.title.y = element_text(
      colour = COL_BLACK,
      margin = margin(r = 10)
    ),
    
    panel.spacing.x = grid::unit(
      1.5,
      "lines"
    ),
    
    panel.spacing.y = grid::unit(
      1.5,
      "lines"
    ),
    
    plot.margin = margin(
      t = 8,
      r = 16,
      b = 16,
      l = 10
    )
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    shape = "none",
    linetype = "none"
  )


save_plot(
  figA2_bias_error,
  file.path(
    fig_supp_dir,
    "04_figA2_mean_bias_by_error_distribution.pdf"
  ),
  width = 12.5,
  height = 8.2
)

utils::write.csv(
  bias_by_error,
  file.path(data_dir, "04_figA2_mean_bias_by_error_data.csv"),
  row.names = FALSE
)


# ==============================================================================
# 13. Supplementary Figure A3: mean bias by predictor distribution
# ==============================================================================

bias_by_x <- estimation_cell %>%
  filter(estimator_id %in% selected_methods_for_sensitivity) %>%
  group_by(
    x_type, outlier_method,
    estimator_id, estimator_label, estimator_order
  ) %>%
  summarise(
    mean_abs_bias = safe_mean(mean_abs_bias),
    .groups = "drop"
  ) %>%
  add_outlier_display() %>%
  mutate(
    x_label = factor(
      unname(x_labels_plot[as.character(x_type)]),
      levels = unname(x_labels_plot[x_order])
    ),
    estimator_id = factor(
      estimator_id,
      levels = estimator_meta$estimator_id
    )
  )


figA3_bias_x <- ggplot(
  bias_by_x,
  aes(
    x = x_label,
    y = mean_abs_bias,
    group = estimator_id,
    colour = estimator_id,
    shape = estimator_id,
    linetype = estimator_id
  )
) +
  geom_line(
    linewidth = 0.55,
    na.rm = TRUE
  ) +
  geom_point(
    size = 2.0,
    stroke = 0.6,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ outlier_label_plot,
    ncol = 2,
    scales = "free_y"
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(
      base = 10,
      sigma = 0.001
    ),
    breaks = compact_pseudolog_breaks,
    labels = compact_scientific_labels,
    expand = expansion(
      mult = c(0.02, 0.08)
    )
  ) +
  scale_colour_manual(
    values = method_colors,
    breaks = MAIN_ESTIMATORS,
    labels = estimator_meta$estimator_label[
      match(
        MAIN_ESTIMATORS,
        estimator_meta$estimator_id
      )
    ],
    drop = TRUE
  ) +
  scale_shape_manual(
    values = method_shapes,
    breaks = MAIN_ESTIMATORS,
    labels = estimator_meta$estimator_label[
      match(
        MAIN_ESTIMATORS,
        estimator_meta$estimator_id
      )
    ],
    drop = TRUE
  ) +
  scale_linetype_manual(
    values = method_linetypes,
    breaks = MAIN_ESTIMATORS,
    labels = estimator_meta$estimator_label[
      match(
        MAIN_ESTIMATORS,
        estimator_meta$estimator_id
      )
    ],
    drop = TRUE
  ) +
  labs(
    x = "Predictor distribution",
    y = "Mean absolute bias",
    colour = "Estimator",
    shape = "Estimator",
    linetype = "Estimator"
  ) +
  theme_paper() +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    shape = "none",
    linetype = "none"
  )


save_plot(
  figA3_bias_x,
  file.path(
    fig_supp_dir,
    "04_figA3_mean_bias_by_x_distribution.pdf"
  ),
  width = 10.6,
  height = 7.1
)

utils::write.csv(
  bias_by_x,
  file.path(data_dir, "04_figA3_mean_bias_by_x_data.csv"),
  row.names = FALSE
)

# ==============================================================================
# 14. Supplementary Figure A4: runtime distributions
# ==============================================================================

runtime_plot_data <- runtime_long %>%
  filter(
    estimator_id %in% MAIN_ESTIMATORS,
    is.finite(runtime_seconds),
    runtime_seconds > 0
  ) %>%
  add_outlier_display() %>%
  mutate(
    estimator_label = factor(
      estimator_label,
      levels = estimator_meta$estimator_label[
        estimator_meta$estimator_id %in% MAIN_ESTIMATORS
      ]
    )
  )


runtime_distribution_summary <- runtime_plot_data %>%
  group_by(
    outlier_method,
    outlier_label_plot,
    estimator_id,
    estimator_label,
    estimator_order
  ) %>%
  summarise(
    n_valid = n(),
    mean_runtime = safe_mean(runtime_seconds),
    sd_runtime = safe_sd(runtime_seconds),
    .groups = "drop"
  )


figA4_runtime <- ggplot(
  runtime_plot_data,
  aes(
    x = runtime_seconds,
    y = estimator_label
  )
) +
  geom_violin(
    fill = COL_VIOLIN,
    colour = COL_GREY_DARK,
    linewidth = 0.25,
    scale = "width",
    trim = TRUE,
    na.rm = TRUE
  ) +
  geom_segment(
    data = runtime_distribution_summary,
    aes(
      x = pmax(
        .Machine$double.eps,
        mean_runtime - sd_runtime
      ),
      xend = mean_runtime + sd_runtime,
      y = estimator_label,
      yend = estimator_label
    ),
    inherit.aes = FALSE,
    colour = COL_BLUE,
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  geom_point(
    data = runtime_distribution_summary,
    aes(
      x = mean_runtime,
      y = estimator_label
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.2,
    stroke = 0.45,
    fill = COL_BLUE,
    colour = COL_BLACK,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ outlier_label_plot,
    ncol = 2
  ) +
  scale_x_log10(
    labels = scales::label_number(
      accuracy = 0.001
    ),
    expand = expansion(
      mult = c(0.02, 0.05)
    )
  ) +
  labs(
    x = "Runtime per Monte Carlo draw (seconds, log scale)",
    y = NULL
  ) +
  theme_distribution() +
  theme(
    legend.position = "none"
  )


save_plot(
  figA4_runtime,
  file.path(
    fig_supp_dir,
    "04_figA4_runtime_distributions.pdf"
  ),
  width = 11.2,
  height = 7.2
)


utils::write.csv(
  runtime_distribution_summary,
  file.path(
    data_dir,
    "04_figA4_runtime_distribution_summary.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 15. Supplementary Figure A5: MIS-oracle bias advantage by DGP
# ==============================================================================

oracle_advantage_heatmap <- estimation_cell %>%
  filter(
    outlier_method != "none",
    estimator_id %in% c("full", "mis_oracle"),
    is.finite(mean_abs_bias),
    mean_abs_bias > 0
  ) %>%
  select(
    x_type,
    error_type,
    outlier_method,
    estimator_id,
    mean_abs_bias
  ) %>%
  pivot_wider(
    names_from = estimator_id,
    values_from = mean_abs_bias
  ) %>%
  filter(
    is.finite(full),
    is.finite(mis_oracle),
    full > 0,
    mis_oracle > 0
  ) %>%
  mutate(
    log2_bias_ratio = log2(
      full / mis_oracle
    ),
    
    predictor_label = factor(
      unname(
        x_labels_table[
          as.character(x_type)
        ]
      ),
      levels = rev(
        unname(
          x_labels_table[x_order]
        )
      )
    ),
    
    error_label = factor(
      unname(
        error_labels_table[
          as.character(error_type)
        ]
      ),
      levels = unname(
        error_labels_table[
          error_order[
            error_order != "gpd"
          ]
        ]
      )
    ),
    
    cell_label = sprintf(
      "%+.1f",
      log2_bias_ratio
    )
  ) %>%
  add_outlier_display()


advantage_limit <- safe_quantile(
  abs(
    oracle_advantage_heatmap$log2_bias_ratio
  ),
  0.95
)

if (
  !is.finite(advantage_limit) ||
  advantage_limit <= 0
) {
  advantage_limit <- 1
}


figA5_oracle_advantage <- ggplot(
  oracle_advantage_heatmap,
  aes(
    x = error_label,
    y = predictor_label,
    fill = log2_bias_ratio
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = cell_label),
    size = 2.7,
    colour = COL_BLACK,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ outlier_label_plot,
    ncol = 3
  ) +
  scale_fill_gradient2(
    low = COL_ORANGE,
    mid = COL_NEUTRAL,
    high = COL_BLUE,
    midpoint = 0,
    limits = c(
      -advantage_limit,
      advantage_limit
    ),
    oob = scales::squish,
    breaks = scales::breaks_pretty(
      n = 5
    ),
    labels = scales::label_number(
      accuracy = 0.1
    ),
    name = "log2(OLS MAB /\nMIS-oracle MAB)"
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  scale_y_discrete(
    drop = FALSE
  ) +
  labs(
    x = "Error distribution",
    y = "Predictor distribution",
    caption = paste(
      "Positive values favour MIS-oracle;",
      "+1 means OLS mean absolute bias is twice",
      "MIS-oracle mean absolute bias."
    )
  ) +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      barwidth = grid::unit(
        8.5,
        "cm"
      ),
      barheight = grid::unit(
        0.55,
        "cm"
      ),
      ticks = TRUE
    )
  ) +
  theme_heatmap() +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    plot.caption = element_text(
      hjust = 0,
      size = 8,
      lineheight = 1.05,
      margin = margin(t = 8)
    )
  )


save_plot(
  figA5_oracle_advantage,
  file.path(
    fig_supp_dir,
    "04_figA5_mis_oracle_bias_advantage_heatmap.pdf"
  ),
  width = 11.5,
  height = 5.8
)


utils::write.csv(
  oracle_advantage_heatmap,
  file.path(
    data_dir,
    "04_figA5_mis_oracle_bias_advantage_data.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 16. Supplementary Figure A6:
#     MIS-oracle versus classical deletion under bad leverage
# ==============================================================================

classical_advantage_heatmap <- estimation_cell %>%
  filter(
    outlier_method == "bad_leverage",
    estimator_id %in% c(
      "cd",
      "dfb",
      "mis_oracle"
    ),
    is.finite(mean_abs_bias),
    mean_abs_bias > 0
  ) %>%
  
  # If future simulation grids contain more than one N or k,
  # give each design cell equal weight within each displayed DGP.
  group_by(
    x_type,
    error_type,
    estimator_id
  ) %>%
  summarise(
    mean_abs_bias =
      safe_mean(mean_abs_bias),
    .groups = "drop"
  ) %>%
  
  pivot_wider(
    names_from = estimator_id,
    values_from = mean_abs_bias
  ) %>%
  
  filter(
    is.finite(cd),
    is.finite(dfb),
    is.finite(mis_oracle),
    cd > 0,
    dfb > 0,
    mis_oracle > 0
  ) %>%
  
  transmute(
    x_type,
    error_type,
    
    `Cook's D` =
      log2(
        cd / mis_oracle
      ),
    
    DFBETAS =
      log2(
        dfb / mis_oracle
      )
  ) %>%
  
  pivot_longer(
    cols = all_of(
      c(
        "Cook's D",
        "DFBETAS"
      )
    ),
    names_to = "benchmark_label",
    values_to = "log2_mae_ratio"
  ) %>%
  
  mutate(
    benchmark_label = factor(
      benchmark_label,
      levels = c(
        "Cook's D",
        "DFBETAS"
      )
    ),
    
    predictor_label = factor(
      unname(
        x_labels_table[
          as.character(x_type)
        ]
      ),
      levels = rev(
        unname(
          x_labels_table[x_order]
        )
      )
    ),
    
    error_label = factor(
      unname(
        error_labels_table[
          as.character(error_type)
        ]
      ),
      levels = unname(
        error_labels_table[
          error_order[
            error_order != "gpd"
          ]
        ]
      )
    ),
    
    cell_label = sprintf(
      "%+.1f",
      log2_mae_ratio
    )
  )


classical_advantage_limit <- safe_quantile(
  abs(
    classical_advantage_heatmap$
      log2_mae_ratio
  ),
  0.95
)


if (
  !is.finite(classical_advantage_limit) ||
  classical_advantage_limit <= 0
) {
  classical_advantage_limit <- 1
}


figA6_classical_advantage <- ggplot(
  classical_advantage_heatmap,
  aes(
    x = error_label,
    y = predictor_label,
    fill = log2_mae_ratio
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  
  geom_text(
    aes(
      label = cell_label
    ),
    size = 2.7,
    colour = COL_BLACK,
    na.rm = TRUE
  ) +
  
  facet_wrap(
    ~ benchmark_label,
    ncol = 2
  ) +
  
  scale_fill_gradient2(
    low = COL_ORANGE,
    mid = COL_NEUTRAL,
    high = COL_BLUE,
    midpoint = 0,
    limits = c(
      -classical_advantage_limit,
      classical_advantage_limit
    ),
    oob = scales::squish,
    breaks = scales::breaks_pretty(
      n = 5
    ),
    labels = scales::label_number(
      accuracy = 0.1
    ),
    name = paste0(
      "log2(classical MAE /\n",
      "MIS-oracle MAE)"
    )
  ) +
  
  scale_x_discrete(
    drop = FALSE
  ) +
  
  scale_y_discrete(
    drop = FALSE
  ) +
  
  labs(
    x = "Error distribution",
    y = "Predictor distribution",
    caption = paste(
      "Bad-leverage contamination only.",
      "Positive values favour MIS-oracle;",
      "+1 means the classical method's MAE is twice",
      "the MIS-oracle MAE."
    )
  ) +
  
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      barwidth = grid::unit(
        8.5,
        "cm"
      ),
      barheight = grid::unit(
        0.55,
        "cm"
      ),
      ticks = TRUE
    )
  ) +
  
  theme_heatmap() +
  
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    
    plot.caption = element_text(
      hjust = 0,
      size = 8,
      lineheight = 1.05,
      margin = margin(t = 8)
    )
  )


save_plot(
  figA6_classical_advantage,
  file.path(
    fig_supp_dir,
    "04_figA6_mis_oracle_vs_classical_bad_leverage_heatmap.pdf"
  ),
  width = 10.8,
  height = 5.2
)


utils::write.csv(
  classical_advantage_heatmap,
  file.path(
    data_dir,
    "04_figA6_mis_oracle_vs_classical_bad_leverage_data.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 19. Supplementary tables
# ==============================================================================

# A1: mean runtime.
tabA1_runtime <- make_metric_wide_table(
  data = runtime_broad %>%
    filter(
      estimator_id %in% MAIN_ESTIMATORS
    ),
  mean_column = "mean_runtime",
  mcse_column = "mcse_runtime",
  formatter = function(x, se) fmt_mean_mcse(x, se, digits = 3L),
  include_clean = TRUE
)

write_tex_table(
  data = tabA1_runtime,
  tex_path = file.path(
    tab_supp_dir,
    "04_tabA1_mean_runtime.tex"
  ),
  caption = paste0(
    "Mean complete runtime in seconds by estimator and contamination mechanism. ",
    "Monte Carlo standard errors are in parentheses."
  ),
  label = "tab:robust-runtime",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrr"
)


# A2: method health and availability.
method_health <- coefficient_long %>%
  filter(
    estimator_id %in% MAIN_ESTIMATORS
  ) %>%
  group_by(
    estimator_id, estimator_label, estimator_order
  ) %>%
  summarise(
    total_rows = n(),
    finite_coefficient_rate = mean(is.finite(coefficient)),
    .groups = "drop"
  ) %>%
  left_join(
    coverage_long %>%
      group_by(estimator_id) %>%
      summarise(
        finite_coverage_rate = mean(is.finite(coverage)),
        .groups = "drop"
      ),
    by = "estimator_id"
  ) %>%
  left_join(
    runtime_long %>%
      group_by(estimator_id) %>%
      summarise(
        finite_runtime_rate = mean(is.finite(runtime_seconds)),
        .groups = "drop"
      ),
    by = "estimator_id"
  ) %>%
  arrange(estimator_order)


tabA2_health <- method_health %>%
  transmute(
    Estimator = estimator_label,
    `Total rows` = total_rows,
    `Finite coefficient` = fmt_pct(
      finite_coefficient_rate,
      digits = 1L
    ),
    `Finite coverage` = fmt_pct(
      finite_coverage_rate,
      digits = 1L
    ),
    `Finite runtime` = fmt_pct(
      finite_runtime_rate,
      digits = 1L
    )
  )


write_tex_table(
  data = tabA2_health,
  tex_path = file.path(
    tab_supp_dir,
    "04_tabA2_method_health.tex"
  ),
  caption = "Method-level availability and numerical health checks.",
  label = "tab:robust-method-health",
  resize_width = TABLE_WIDTH_COMPACT,
  align = "lrrrr"
)

# A3: bad-leverage DGP-level estimation and detection detail.

bad_leverage_mae <- estimation_cell %>%
  filter(
    outlier_method == "bad_leverage",
    estimator_id %in% MAIN_ESTIMATORS,
    is.finite(mean_abs_bias)
  ) %>%
  group_by(
    x_type,
    error_type,
    estimator_id
  ) %>%
  summarise(
    MAE =
      safe_mean(mean_abs_bias),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = estimator_id,
    values_from = MAE
  )


bad_leverage_precision <- detection_quality_cell %>%
  filter(
    outlier_method == "bad_leverage",
    estimator_id %in% MAIN_SELECTION_METHODS,
    is.finite(mean_precision)
  ) %>%
  group_by(
    x_type,
    error_type,
    estimator_id
  ) %>%
  summarise(
    precision =
      safe_mean(mean_precision),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = estimator_id,
    values_from = precision,
    names_prefix = "precision_"
  )


bad_leverage_detail <- bad_leverage_mae %>%
  left_join(
    bad_leverage_precision,
    by = c(
      "x_type",
      "error_type"
    )
  ) %>%
  mutate(
    predictor_label = factor(
      unname(
        x_labels_table[
          as.character(x_type)
        ]
      ),
      levels = unname(
        x_labels_table[x_order]
      )
    ),
    
    error_label = factor(
      unname(
        error_labels_table[
          as.character(error_type)
        ]
      ),
      levels = unname(
        error_labels_table[
          error_order[
            error_order != "gpd"
          ]
        ]
      )
    )
  ) %>%
  arrange(
    predictor_label,
    error_label
  )


tabA3_bad_leverage <- bad_leverage_detail %>%
  transmute(
    Predictor =
      as.character(predictor_label),
    
    Error =
      as.character(error_label),
    
    `OLS MAE` =
      fmt_num(
        full,
        digits = 5L
      ),
    
    `Cook MAE` =
      fmt_num(
        cd,
        digits = 5L
      ),
    
    `DFBETAS MAE` =
      fmt_num(
        dfb,
        digits = 5L
      ),
    
    `MIS-oracle MAE` =
      fmt_num(
        mis_oracle,
        digits = 5L
      ),
    
    `MM MAE` =
      fmt_num(
        mm,
        digits = 5L
      ),
    
    `LTS MAE` =
      fmt_num(
        lts,
        digits = 5L
      ),
    
    `Cook precision` =
      fmt_pct(
        precision_cd,
        digits = 1L
      ),
    
    `DFBETAS precision` =
      fmt_pct(
        precision_dfb,
        digits = 1L
      ),
    
    `MIS-oracle precision` =
      fmt_pct(
        precision_mis_oracle,
        digits = 1L
      )
  )


write_tex_table(
  data = tabA3_bad_leverage,
  tex_path = file.path(
    tab_supp_dir,
    "04_tabA3_bad_leverage_dgp_detail.tex"
  ),
  caption = paste0(
    "Bad-leverage performance by predictor and error distribution. ",
    "All cells use N = 5000, true coalition size k = 50, and a ",
    "1 percent injected coalition. MAE is mean absolute coefficient error. ",
    "Precision is the fraction of selected observations belonging to the ",
    "injected coalition. Because MIS with oracle k selects exactly 50 ",
    "observations when the true coalition also contains 50 observations, ",
    "its precision equals its recall. The legacy GPD condition is excluded ",
    "from this finite-moment comparison."
  ),
  label = "tab:robust-bad-leverage-dgp-detail",
  resize_width = TABLE_WIDTH_WIDE,
  align = "llrrrrrrrrr",
  font_command = "\\scriptsize"
)


utils::write.csv(
  bad_leverage_detail,
  file.path(
    data_dir,
    "04_tabA3_bad_leverage_dgp_detail.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 20. Additional machine-readable full results
# ==============================================================================

# Complete estimator-cell table.
complete_cell_table <- estimation_cell %>%
  arrange(
    outlier_method,
    n_obs,
    contam_prop,
    x_type,
    error_type,
    estimator_order
  )

utils::write.csv(
  complete_cell_table,
  file.path(
    data_dir,
    "04_complete_estimation_results_by_cell.csv"
  ),
  row.names = FALSE
)


# Complete selected-k and overlap table.
complete_detection_cell_table <- full_join(
  selection_cell,
  overlap_cell,
  by = c(
    "n_obs",
    "design_k",
    "contam_prop",
    "x_type",
    "error_type",
    "outlier_method",
    "estimator_id",
    "estimator_label",
    "estimator_order",
    "selection_order"
  )
) %>%
  arrange(
    outlier_method,
    n_obs,
    contam_prop,
    x_type,
    error_type,
    selection_order
  )

utils::write.csv(
  complete_detection_cell_table,
  file.path(
    data_dir,
    "04_complete_selection_overlap_by_cell.csv"
  ),
  row.names = FALSE
)


# Save all publication summaries in one reusable RDS.
publication_summaries <- list(
  estimation_cell = estimation_cell,
  estimation_broad = estimation_broad,
  
  overlap_cell = overlap_cell,
  overlap_broad = overlap_broad,
  
  runtime_cell = runtime_cell,
  runtime_broad = runtime_broad,
  
  method_health = method_health,
  
  figure1_summary = coef_distribution_summary,
  
  figure2_data = tail_mis_oracle_data,
  figure2_mis_oracle_data = tail_mis_oracle_data,
  figure2_full_tail_summary = tail_cell_summary,
  figure2_problem_dgp_cells = problem_dgp_cells,
  
  figure3_summary = detection_quality_distribution_summary,
  
  figure3_cell = detection_quality_cell,
  
  figure3_broad = detection_quality_broad,
  
  figureA4_runtime_summary = runtime_distribution_summary,
  
  figureA5_oracle_advantage_data = oracle_advantage_heatmap,
  
  figureA6_classical_advantage_data = classical_advantage_heatmap,
  
  tableA3_bad_leverage_detail = bad_leverage_detail
)

saveRDS( publication_summaries,
  file.path(data_dir, "04_publication_summaries.rds")
)


# ==============================================================================
# 21. Diagnostics, clipping audit, and manifest
# ==============================================================================

n_coef_finite_displayed <- nrow(
  coef_finite
)

n_coef_nonfinite_excluded <- coefficient_long %>%
  filter(
    estimator_id %in% MAIN_ESTIMATORS
  ) %>%
  summarise(
    n = sum(
      !is.finite(signed_error)
    )
  ) %>%
  pull(n)

input_audit <- data.frame(
  item = c(
    "Rows in primary RDS",
    "Unique x distributions",
    "Unique error distributions",
    "Unique contamination mechanisms",
    "Unique Monte Carlo iteration IDs",
    "Maximum recorded-versus-recomputed bias difference",
    "Finite iteration-level coefficient errors displayed in Figure 1",
    "Non-finite iteration-level coefficient errors excluded from Figure 1",
    "Optional summary RDS present",
    "Optional bias-summary RDS present"
  ),
  value = c(
    as.character(nrow(sim)),
    as.character(length(unique(sim$x_type))),
    as.character(length(unique(sim$error_type))),
    as.character(length(unique(sim$outlier_method))),
    as.character(length(unique(sim$iter))),
    as.character(max_bias_difference),
    as.character(n_coef_finite_displayed),
    as.character(n_coef_nonfinite_excluded),
    ifelse(file.exists(input_summary_optional), "yes", "no"),
    ifelse(file.exists(input_bias_optional), "yes", "no")
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  input_audit,
  file.path(diag_dir, "04_input_and_plot_audit.csv"),
  row.names = FALSE
)


n_detection_quality_outside <- sum(
  detection_quality_plot_data$metric_value < 0 |
    detection_quality_plot_data$metric_value > 1,
  na.rm = TRUE
)


clipping_audit <- data.frame(
  figure = c(
    "04_fig1_coefficient_error_distributions",
    "04_fig3_detection_quality_distributions"
  ),
  
  lower_limit = c(
    NA_real_,
    0
  ),
  
  upper_limit = c(
    NA_real_,
    1
  ),
  
  observations_outside_display = c(
    0L,
    n_detection_quality_outside
  ),
  
  note = c(
    "Figure 1 uses a pseudo-log transformation without coordinate clipping.",
    "Detection recall and precision are displayed on their complete [0, 1] range."
  ),
  
  stringsAsFactors = FALSE
)

utils::write.csv(
  clipping_audit,
  file.path(diag_dir, "04_clipping_audit.csv"),
  row.names = FALSE
)

cat(
  "Main coefficient figure uses a signed pseudo-logarithmic axis ",
  "with no clipping of finite coefficient errors.\n",
  sep = ""
)

cat(sprintf(
  "Finite coefficient errors displayed in Figure 1: %d\n",
  n_coef_finite_displayed
))

cat(sprintf(
  "Non-finite coefficient errors excluded from Figure 1: %d\n",
  n_coef_nonfinite_excluded
))

cat(
  "Arithmetic means remain the primary summaries; ",
  "medians in Figure 1 are secondary markers of distributional skewness.\n",
  sep = ""
)