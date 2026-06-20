# ==============================================================================
# File: scripts/84_output_robust.R
# Purpose:
#   Generate publication-ready figures and tables for Script 04:
#   comparison of OLS, classical diagnostic deletion, MIS variants,
#   selection-adjusted permutation MIS (MIS-SAP), MM regression, and LTS.
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
#   \usepackage{float}
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
#   - Full empirical distributions are shown with violin densities.
#   - Arithmetic means are marked explicitly; no boxplots or median lines.
#   - Labels, legends, facets, and margins are sized for formal papers.
#   - No red-green palette is used.
#   - Diverging heatmaps use orange-white-blue:
#       orange = MIS-SAP worse, blue = MIS-SAP better.
#
# Table design rules:
#   - Tables are complete \begin{table}[H] environments.
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

# Main coefficient figure displays the central 99% of finite estimates.
# A full-range coefficient-distribution figure is also saved in the supplement.
COEF_MAIN_QUANTILES <- c(0.005, 0.995)

# Main selected-k figure displays the central 99.5% of finite selected-k values.
# A full-range data file remains available.
K_MAIN_UPPER_QUANTILE <- 0.995

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
    placement = "H",
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
    "iter", "x_type", "error_type", "outlier_method", "set_size"
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

required_columns <- c(
  "iter", "x_type", "error_type", "outlier_method", "set_size",
  
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
# 4. Reshape iteration-level data
# ==============================================================================

coefficient_long <- reshape_mapped_columns(
  data = sim,
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
  data = sim,
  mapping = bias_mapping,
  value_name = "absolute_bias_recorded",
  metadata = estimator_meta,
  metadata_key = "estimator_id"
)

coverage_long <- reshape_mapped_columns(
  data = sim,
  mapping = coverage_mapping,
  value_name = "coverage",
  metadata = estimator_meta,
  metadata_key = "estimator_id"
)

runtime_long <- reshape_mapped_columns(
  data = sim,
  mapping = runtime_mapping,
  value_name = "runtime_seconds",
  metadata = estimator_meta,
  metadata_key = "estimator_id"
)

selection_long <- reshape_mapped_columns(
  data = sim,
  mapping = k_mapping,
  value_name = "selected_k",
  metadata = selection_meta,
  metadata_key = "estimator_id"
)

overlap_long <- reshape_mapped_columns(
  data = sim,
  mapping = overlap_mapping,
  value_name = "overlap",
  metadata = selection_meta,
  metadata_key = "estimator_id"
) %>%
  filter(outlier_method != "none")


# Validate that recorded absolute bias agrees with the coefficient columns.
bias_validation <- coefficient_long %>%
  select(
    iter, x_type, error_type, outlier_method,
    estimator_id, coefficient
  ) %>%
  left_join(
    bias_long %>%
      select(
        iter, x_type, error_type, outlier_method,
        estimator_id, absolute_bias_recorded
      ),
    by = c(
      "iter", "x_type", "error_type",
      "outlier_method", "estimator_id"
    )
  ) %>%
  mutate(
    absolute_bias_recomputed = abs(coefficient - TRUE_BETA),
    difference = absolute_bias_recorded - absolute_bias_recomputed
  )

max_bias_difference <- if (
  any(is.finite(bias_validation$difference))
) {
  max(abs(bias_validation$difference), na.rm = TRUE)
} else {
  NA_real_
}

if (is.finite(max_bias_difference) && max_bias_difference > 1e-8) {
  warning(
    "Recorded bias columns differ from abs(coefficient - TRUE_BETA). ",
    "Maximum absolute difference: ",
    format(max_bias_difference, scientific = TRUE),
    ". Publication summaries use coefficient-based recomputation."
  )
}


# ==============================================================================
# 5. Cell-level and equal-cell-weight broad summaries
# ==============================================================================

estimation_cell <- coefficient_long %>%
  group_by(
    x_type, error_type, outlier_method,
    estimator_id, estimator_label, estimator_short,
    estimator_family, estimator_order
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
    x_type, error_type, outlier_method,
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
      "x_type", "error_type", "outlier_method", "estimator_id"
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
    x_type, error_type, outlier_method,
    estimator_id, estimator_label,
    estimator_order, selection_order
  ) %>%
  summarise(
    n_valid_k = sum(is.finite(selected_k)),
    mean_selected_k = safe_mean(selected_k),
    sd_selected_k = safe_sd(selected_k),
    mcse_selected_k = safe_mcse(selected_k),
    .groups = "drop"
  )


selection_broad <- selection_cell %>%
  group_by(
    outlier_method,
    estimator_id, estimator_label,
    estimator_order, selection_order
  ) %>%
  summarise(
    n_design_cells = sum(is.finite(mean_selected_k)),
    mean_selected_k = safe_mean(mean_selected_k),
    mcse_selected_k = combined_cell_mcse(
      mcse_selected_k,
      mean_selected_k
    ),
    .groups = "drop"
  )


overlap_cell <- overlap_long %>%
  group_by(
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


sap_cell <- sim %>%
  group_by(x_type, error_type, outlier_method) %>%
  summarise(
    n_iter = n(),
    
    detection_rate = safe_mean(as.numeric(k_peel_sap > 0L)),
    mcse_detection = safe_prop_mcse(as.numeric(k_peel_sap > 0L)),
    
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
# 6. Main Figure 1: coefficient-estimate distributions
# ==============================================================================

coef_plot_data <- coefficient_long %>%
  filter(is.finite(coefficient)) %>%
  add_outlier_display() %>%
  mutate(
    estimator_label = factor(
      estimator_label,
      levels = estimator_meta$estimator_label
    )
  )


coef_distribution_summary <- coefficient_long %>%
  group_by(
    outlier_method,
    estimator_id, estimator_label, estimator_order
  ) %>%
  summarise(
    n_valid = sum(is.finite(coefficient)),
    mean_coefficient = safe_mean(coefficient),
    sd_coefficient = safe_sd(coefficient),
    .groups = "drop"
  ) %>%
  add_outlier_display() %>%
  mutate(
    estimator_label = factor(
      estimator_label,
      levels = estimator_meta$estimator_label
    )
  )


coef_limits <- c(
  safe_quantile(
    coef_plot_data$coefficient,
    COEF_MAIN_QUANTILES[[1L]]
  ),
  safe_quantile(
    coef_plot_data$coefficient,
    COEF_MAIN_QUANTILES[[2L]]
  )
)

coef_limits[[1L]] <- min(coef_limits[[1L]], TRUE_BETA, na.rm = TRUE)
coef_limits[[2L]] <- max(coef_limits[[2L]], TRUE_BETA, na.rm = TRUE)

coef_padding <- 0.04 * diff(coef_limits)

if (!is.finite(coef_padding) || coef_padding <= 0) {
  coef_padding <- 0.1
}

coef_limits <- coef_limits + c(-coef_padding, coef_padding)

n_coef_outside_main <- sum(
  coef_plot_data$coefficient < coef_limits[[1L]] |
    coef_plot_data$coefficient > coef_limits[[2L]],
  na.rm = TRUE
)


fig1_coef <- ggplot(
  coef_plot_data,
  aes(x = coefficient, y = estimator_label)
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
    data = coef_distribution_summary,
    aes(
      x = mean_coefficient - sd_coefficient,
      xend = mean_coefficient + sd_coefficient,
      y = estimator_label,
      yend = estimator_label
    ),
    inherit.aes = FALSE,
    colour = COL_BLUE,
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  geom_point(
    data = coef_distribution_summary,
    aes(
      x = mean_coefficient,
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
    xintercept = TRUE_BETA,
    colour = COL_ORANGE,
    linetype = "dashed",
    linewidth = 0.75
  ) +
  facet_wrap(
    ~ outlier_label_plot,
    ncol = 2,
    scales = "fixed"
  ) +
  coord_cartesian(
    xlim = coef_limits,
    clip = "on"
  ) +
  labs(
    x = expression(hat(beta)),
    y = NULL
  ) +
  theme_distribution() +
  theme(
    legend.position = "none"
  )


save_plot(
  fig1_coef,
  file.path(fig_main_dir, "04_fig1_coefficient_distributions.pdf"),
  width = 10.4,
  height = 7.6
)

utils::write.csv(
  coef_distribution_summary,
  file.path(data_dir, "04_fig1_coefficient_distribution_summary.csv"),
  row.names = FALSE
)


# Full-range version for supplementary material.
figA0_coef_full <- fig1_coef +
  coord_cartesian(clip = "on")

save_plot(
  figA0_coef_full,
  file.path(
    fig_supp_dir,
    "04_figA0_coefficient_distributions_full_range.pdf"
  ),
  width = 10.4,
  height = 7.6
)


# ==============================================================================
# 7. Main Figure 2: mean absolute bias versus coverage
# ==============================================================================

performance_plot_data <- estimation_broad %>%
  filter(
    is.finite(mean_abs_bias),
    is.finite(coverage_rate)
  ) %>%
  add_outlier_display() %>%
  mutate(
    estimator_id = factor(
      estimator_id,
      levels = estimator_meta$estimator_id
    )
  )


fig2_tradeoff <- ggplot(
  performance_plot_data,
  aes(
    x = mean_abs_bias,
    y = coverage_rate,
    colour = estimator_id,
    shape = estimator_id
  )
) +
  geom_hline(
    yintercept = 0.95,
    colour = COL_GREY_DARK,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_point(
    size = 2.8,
    stroke = 0.7
  ) +
  geom_text(
    aes(label = estimator_short),
    nudge_y = 0.025,
    size = 2.55,
    show.legend = FALSE,
    check_overlap = FALSE
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
    labels = scales::label_number(accuracy = 0.001)
  ) +
  scale_y_continuous(
    limits = c(0, 1.06),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1)
  ) +
  scale_colour_manual(
    values = method_colors,
    breaks = estimator_meta$estimator_id,
    labels = estimator_meta$estimator_label,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = method_shapes,
    breaks = estimator_meta$estimator_id,
    labels = estimator_meta$estimator_label,
    drop = FALSE
  ) +
  labs(
    x = "Mean absolute bias",
    y = "Empirical 95% CI coverage",
    colour = "Estimator",
    shape = "Estimator"
  ) +
  theme_paper() +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 14, 10, 10)
  )


save_plot(
  fig2_tradeoff,
  file.path(fig_main_dir, "04_fig2_bias_coverage_tradeoff.pdf"),
  width = 10.4,
  height = 7.4
)

utils::write.csv(
  performance_plot_data,
  file.path(data_dir, "04_fig2_bias_coverage_tradeoff_data.csv"),
  row.names = FALSE
)


# ==============================================================================
# 8. Main Figure 3: selected-k distributions
# ==============================================================================

selection_plot_data <- selection_long %>%
  filter(is.finite(selected_k)) %>%
  add_outlier_display() %>%
  mutate(
    estimator_label = factor(
      estimator_label,
      levels = selection_meta$estimator_label[
        order(selection_meta$selection_order)
      ]
    )
  )


selection_distribution_summary <- selection_long %>%
  group_by(
    outlier_method,
    estimator_id, estimator_label, selection_order
  ) %>%
  summarise(
    n_valid = sum(is.finite(selected_k)),
    mean_selected_k = safe_mean(selected_k),
    sd_selected_k = safe_sd(selected_k),
    .groups = "drop"
  ) %>%
  add_outlier_display() %>%
  mutate(
    estimator_label = factor(
      estimator_label,
      levels = selection_meta$estimator_label[
        order(selection_meta$selection_order)
      ]
    )
  )


true_k_reference <- sim %>%
  group_by(outlier_method) %>%
  summarise(
    true_k = safe_mean(set_size),
    .groups = "drop"
  ) %>%
  add_outlier_display()


k_upper <- safe_quantile(
  selection_plot_data$selected_k,
  K_MAIN_UPPER_QUANTILE
)

k_upper <- max(
  k_upper,
  max(true_k_reference$true_k, na.rm = TRUE),
  na.rm = TRUE
)

if (!is.finite(k_upper) || k_upper <= 0) {
  k_upper <- 1
}

k_upper <- k_upper * 1.05

n_k_outside_main <- sum(
  selection_plot_data$selected_k > k_upper,
  na.rm = TRUE
)


fig3_k <- ggplot(
  selection_plot_data,
  aes(x = selected_k, y = estimator_label)
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
    data = selection_distribution_summary,
    aes(
      x = pmax(0, mean_selected_k - sd_selected_k),
      xend = mean_selected_k + sd_selected_k,
      y = estimator_label,
      yend = estimator_label
    ),
    inherit.aes = FALSE,
    colour = COL_BLUE,
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  geom_point(
    data = selection_distribution_summary,
    aes(
      x = mean_selected_k,
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
    data = true_k_reference,
    aes(xintercept = true_k),
    colour = COL_ORANGE,
    linetype = "dashed",
    linewidth = 0.75
  ) +
  facet_wrap(
    ~ outlier_label_plot,
    ncol = 2,
    scales = "fixed"
  ) +
  coord_cartesian(
    xlim = c(0, k_upper),
    clip = "on"
  ) +
  labs(
    x = "Selected number of observations",
    y = NULL
  ) +
  theme_distribution() +
  theme(
    legend.position = "none"
  )


save_plot(
  fig3_k,
  file.path(fig_main_dir, "04_fig3_selected_k_distributions.pdf"),
  width = 10.4,
  height = 7.1
)

utils::write.csv(
  selection_distribution_summary,
  file.path(data_dir, "04_fig3_selected_k_distribution_summary.csv"),
  row.names = FALSE
)


# ==============================================================================
# 9. Main Figure 4: detection-overlap distributions
# ==============================================================================

overlap_plot_data <- overlap_long %>%
  filter(is.finite(overlap)) %>%
  add_outlier_display() %>%
  mutate(
    estimator_label = factor(
      estimator_label,
      levels = selection_meta$estimator_label[
        order(selection_meta$selection_order)
      ]
    )
  )


overlap_distribution_summary <- overlap_long %>%
  group_by(
    outlier_method,
    estimator_id, estimator_label, selection_order
  ) %>%
  summarise(
    n_valid = sum(is.finite(overlap)),
    mean_overlap = safe_mean(overlap),
    sd_overlap = safe_sd(overlap),
    .groups = "drop"
  ) %>%
  add_outlier_display() %>%
  mutate(
    estimator_label = factor(
      estimator_label,
      levels = selection_meta$estimator_label[
        order(selection_meta$selection_order)
      ]
    )
  )


fig4_overlap <- ggplot(
  overlap_plot_data,
  aes(x = overlap, y = estimator_label)
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
    data = overlap_distribution_summary,
    aes(
      x = pmax(0, mean_overlap - sd_overlap),
      xend = pmin(1, mean_overlap + sd_overlap),
      y = estimator_label,
      yend = estimator_label
    ),
    inherit.aes = FALSE,
    colour = COL_BLUE,
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  geom_point(
    data = overlap_distribution_summary,
    aes(
      x = mean_overlap,
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
  facet_wrap(
    ~ outlier_label_plot,
    ncol = 3
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1)
  ) +
  labs(
    x = "Fraction of injected observations recovered",
    y = NULL
  ) +
  theme_distribution() +
  theme(
    legend.position = "none"
  )


save_plot(
  fig4_overlap,
  file.path(fig_main_dir, "04_fig4_detection_overlap_distributions.pdf"),
  width = 11.2,
  height = 5.1
)

utils::write.csv(
  overlap_distribution_summary,
  file.path(data_dir, "04_fig4_overlap_distribution_summary.csv"),
  row.names = FALSE
)


# ==============================================================================
# 10. Main tables
# ==============================================================================

# Table 1: arithmetic mean coefficients.
tab1_mean_coef <- make_metric_wide_table(
  data = estimation_broad,
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
    "mechanism. Monte Carlo standard errors are in parentheses. ",
    "The true coefficient is beta0 = 1."
  ),
  label = "tab:robust-mean-coefficients",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrr"
)


# Table 2a: mean absolute bias.
tab2a_bias <- make_metric_wide_table(
  data = estimation_broad,
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
    "Monte Carlo standard errors are in parentheses."
  ),
  label = "tab:robust-mean-absolute-bias",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrr"
)


# Table 2b: RMSE.
tab2b_rmse <- make_metric_wide_table(
  data = estimation_broad,
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
    "Delta-method Monte Carlo standard errors are in parentheses."
  ),
  label = "tab:robust-rmse",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrr"
)


# Table 2c: empirical coverage.
tab2c_coverage <- make_metric_wide_table(
  data = estimation_broad,
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
    "Monte Carlo standard errors in percentage points are in parentheses."
  ),
  label = "tab:robust-coverage",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrr",
  escape_cells = FALSE
)


# Table 3a: selected k.
tab3a_selected_k <- make_metric_wide_table(
  data = selection_broad,
  mean_column = "mean_selected_k",
  mcse_column = "mcse_selected_k",
  formatter = function(x, se) fmt_mean_mcse(x, se, digits = 2L),
  include_clean = TRUE
)

write_tex_table(
  data = tab3a_selected_k,
  tex_path = file.path(
    tab_main_dir,
    "04_tab3a_mean_selected_k.tex"
  ),
  caption = paste0(
    "Mean number of observations selected or removed by each detection method. ",
    "Monte Carlo standard errors are in parentheses."
  ),
  label = "tab:robust-selected-k",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrr"
)


# Table 3b: detection overlap.
tab3b_overlap <- make_metric_wide_table(
  data = overlap_broad,
  mean_column = "mean_overlap",
  mcse_column = "mcse_overlap",
  formatter = function(x, se) fmt_pct_mcse(x, se, digits = 1L),
  include_clean = FALSE
)

write_tex_table(
  data = tab3b_overlap,
  tex_path = file.path(
    tab_main_dir,
    "04_tab3b_mean_detection_overlap.tex"
  ),
  caption = paste0(
    "Mean fraction of injected observations recovered by each detection method. ",
    "Monte Carlo standard errors in percentage points are in parentheses."
  ),
  label = "tab:robust-detection-overlap",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrr",
  escape_cells = FALSE
)


# Table 4: SAP process and inferential diagnostics.
tab4_sap <- sap_broad %>%
  add_outlier_display() %>%
  transmute(
    Scenario = as.character(outlier_label_table),
    `Detection rate` = fmt_pct_mcse(
      detection_rate,
      mcse_detection,
      digits = 1L
    ),
    `Mean selected k` = fmt_mean_mcse(
      mean_selected_k,
      mcse_selected_k,
      digits = 2L
    ),
    `Exact-k rate` = fmt_pct_mcse(
      exact_k_rate,
      mcse_exact_k,
      digits = 1L
    ),
    `Mean overlap` = fmt_pct_mcse(
      mean_overlap,
      mcse_overlap,
      digits = 1L
    ),
    `Mean final p` = fmt_mean_mcse(
      mean_final_p,
      mcse_final_p,
      digits = 3L
    ),
    `Error rate` = fmt_pct_mcse(
      error_rate,
      mcse_error,
      digits = 1L
    ),
    `Mean abs. bias` = fmt_mean_mcse(
      mean_abs_bias_sap,
      mcse_abs_bias_sap,
      digits = 3L
    ),
    Coverage = fmt_pct_mcse(
      coverage_sap,
      mcse_coverage_sap,
      digits = 1L
    )
  )

write_tex_table(
  data = tab4_sap,
  tex_path = file.path(
    tab_main_dir,
    "04_tab4_sap_diagnostics.tex"
  ),
  caption = paste0(
    "Selection-adjusted permutation MIS diagnostics by contamination mechanism. ",
    "Entries are arithmetic means with Monte Carlo standard errors in parentheses."
  ),
  label = "tab:robust-sap-diagnostics",
  resize_width = TABLE_WIDTH_WIDE,
  align = "lrrrrrrrr",
  escape_cells = FALSE,
  font_command = "\\small"
)


# ==============================================================================
# 11. Supplementary Figure A1: absolute-bias distributions
# ==============================================================================

bias_plot_data <- coefficient_long %>%
  mutate(absolute_bias = abs(coefficient - TRUE_BETA)) %>%
  filter(is.finite(absolute_bias)) %>%
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
    labels = scales::label_number(accuracy = 0.001)
  ) +
  labs(
    x = "Absolute coefficient bias",
    y = NULL
  ) +
  theme_distribution() +
  theme(
    legend.position = "none"
  )


save_plot(
  figA1_bias,
  file.path(
    fig_supp_dir,
    "04_figA1_absolute_bias_distributions.pdf"
  ),
  width = 10.4,
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

selected_methods_for_sensitivity <- c(
  "full",
  "dfb",
  "mis_alpha",
  "mis_peel",
  "mis_sap",
  "mis_oracle",
  "mm",
  "lts"
)


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
  scale_y_continuous(
    trans = scales::pseudo_log_trans(
      base = 10,
      sigma = 0.001
    ),
    labels = scales::label_number(accuracy = 0.001)
  ) +
  scale_colour_manual(
    values = method_colors,
    breaks = estimator_meta$estimator_id,
    labels = estimator_meta$estimator_label,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = method_shapes,
    breaks = estimator_meta$estimator_id,
    labels = estimator_meta$estimator_label,
    drop = FALSE
  ) +
  scale_linetype_manual(
    values = method_linetypes,
    breaks = estimator_meta$estimator_id,
    labels = estimator_meta$estimator_label,
    drop = FALSE
  ) +
  labs(
    x = "Error distribution",
    y = "Mean absolute bias",
    colour = "Estimator",
    shape = "Estimator",
    linetype = "Estimator"
  ) +
  theme_paper() +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1
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
  width = 11.0,
  height = 7.5
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
    labels = scales::label_number(accuracy = 0.001)
  ) +
  scale_colour_manual(
    values = method_colors,
    breaks = estimator_meta$estimator_id,
    labels = estimator_meta$estimator_label,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = method_shapes,
    breaks = estimator_meta$estimator_id,
    labels = estimator_meta$estimator_label,
    drop = FALSE
  ) +
  scale_linetype_manual(
    values = method_linetypes,
    breaks = estimator_meta$estimator_id,
    labels = estimator_meta$estimator_label,
    drop = FALSE
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
    is.finite(runtime_seconds),
    runtime_seconds >= 0
  ) %>%
  mutate(
    estimator_label = factor(
      estimator_label,
      levels = estimator_meta$estimator_label
    )
  )


runtime_distribution_summary <- runtime_plot_data %>%
  group_by(
    estimator_id, estimator_label, estimator_order
  ) %>%
  summarise(
    n_valid = n(),
    mean_runtime = safe_mean(runtime_seconds),
    sd_runtime = safe_sd(runtime_seconds),
    .groups = "drop"
  )


figA4_runtime <- ggplot(
  runtime_plot_data,
  aes(x = runtime_seconds, y = estimator_label)
) +
  geom_violin(
    fill = COL_VIOLIN,
    colour = COL_GREY_DARK,
    linewidth = 0.25,
    scale = "width",
    trim = TRUE
  ) +
  geom_segment(
    data = runtime_distribution_summary,
    aes(
      x = pmax(0, mean_runtime - sd_runtime),
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
  scale_x_continuous(
    trans = scales::pseudo_log_trans(
      base = 10,
      sigma = 1e-4
    ),
    labels = scales::label_number(accuracy = 0.001)
  ) +
  labs(
    x = "Runtime per Monte Carlo draw (seconds)",
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
  width = 8.3,
  height = 5.7
)

utils::write.csv(
  runtime_distribution_summary,
  file.path(data_dir, "04_figA4_runtime_distribution_summary.csv"),
  row.names = FALSE
)


# ==============================================================================
# 15. Supplementary Figure A5: SAP p-value distributions
# ==============================================================================

sap_pvalue_long <- sim %>%
  select(
    iter, x_type, error_type, outlier_method,
    peel_sap_final_p, peel_sap_min_p
  ) %>%
  pivot_longer(
    cols = c(peel_sap_final_p, peel_sap_min_p),
    names_to = "p_metric",
    values_to = "p_value"
  ) %>%
  mutate(
    p_metric_label = factor(
      p_metric,
      levels = c("peel_sap_final_p", "peel_sap_min_p"),
      labels = c("Final global p-value", "Minimum global p-value")
    )
  ) %>%
  filter(is.finite(p_value)) %>%
  add_outlier_display()


sap_pvalue_summary <- sap_pvalue_long %>%
  group_by(
    outlier_method,
    outlier_label_plot,
    p_metric,
    p_metric_label
  ) %>%
  summarise(
    n_valid = n(),
    mean_p = safe_mean(p_value),
    sd_p = safe_sd(p_value),
    .groups = "drop"
  )


figA5_pvalues <- ggplot(
  sap_pvalue_long,
  aes(x = p_value, y = p_metric_label)
) +
  geom_violin(
    fill = COL_VIOLIN,
    colour = COL_GREY_DARK,
    linewidth = 0.25,
    scale = "width",
    trim = TRUE
  ) +
  geom_segment(
    data = sap_pvalue_summary,
    aes(
      x = pmax(0, mean_p - sd_p),
      xend = pmin(1, mean_p + sd_p),
      y = p_metric_label,
      yend = p_metric_label
    ),
    inherit.aes = FALSE,
    colour = COL_BLUE,
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  geom_point(
    data = sap_pvalue_summary,
    aes(
      x = mean_p,
      y = p_metric_label
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
    xintercept = 0.05,
    colour = COL_ORANGE,
    linetype = "dashed",
    linewidth = 0.75
  ) +
  facet_wrap(
    ~ outlier_label_plot,
    ncol = 2
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(
    x = "SAP p-value",
    y = NULL
  ) +
  theme_distribution() +
  theme(
    legend.position = "none"
  )


save_plot(
  figA5_pvalues,
  file.path(
    fig_supp_dir,
    "04_figA5_sap_pvalue_distributions.pdf"
  ),
  width = 10.2,
  height = 5.8
)

utils::write.csv(
  sap_pvalue_summary,
  file.path(data_dir, "04_figA5_sap_pvalue_summary.csv"),
  row.names = FALSE
)


# ==============================================================================
# 16. Supplementary Figure A6: SAP peak-excess distributions
# ==============================================================================

sap_excess_data <- sim %>%
  filter(is.finite(peel_sap_peak_excess)) %>%
  add_outlier_display()


sap_excess_summary <- sap_excess_data %>%
  group_by(outlier_method, outlier_label_plot) %>%
  summarise(
    n_valid = n(),
    mean_peak_excess = safe_mean(peel_sap_peak_excess),
    sd_peak_excess = safe_sd(peel_sap_peak_excess),
    .groups = "drop"
  )


figA6_excess <- ggplot(
  sap_excess_data,
  aes(
    x = peel_sap_peak_excess,
    y = outlier_label_plot
  )
) +
  geom_violin(
    fill = COL_VIOLIN,
    colour = COL_GREY_DARK,
    linewidth = 0.25,
    scale = "width",
    trim = TRUE
  ) +
  geom_segment(
    data = sap_excess_summary,
    aes(
      x = pmax(0, mean_peak_excess - sd_peak_excess),
      xend = mean_peak_excess + sd_peak_excess,
      y = outlier_label_plot,
      yend = outlier_label_plot
    ),
    inherit.aes = FALSE,
    colour = COL_BLUE,
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  geom_point(
    data = sap_excess_summary,
    aes(
      x = mean_peak_excess,
      y = outlier_label_plot
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.2,
    stroke = 0.45,
    fill = COL_BLUE,
    colour = COL_BLACK,
    na.rm = TRUE
  ) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(
      base = 10,
      sigma = 0.01
    ),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  labs(
    x = "SAP peak excess ratio",
    y = NULL
  ) +
  theme_distribution() +
  theme(
    legend.position = "none"
  )


save_plot(
  figA6_excess,
  file.path(
    fig_supp_dir,
    "04_figA6_sap_peak_excess_distributions.pdf"
  ),
  width = 8.3,
  height = 4.5
)

utils::write.csv(
  sap_excess_summary,
  file.path(data_dir, "04_figA6_sap_peak_excess_summary.csv"),
  row.names = FALSE
)


# ==============================================================================
# 17. Supplementary Figure A7: SAP stopping reasons
# ==============================================================================

sap_stop_data <- sim %>%
  mutate(
    peel_sap_stop = ifelse(
      is.na(peel_sap_stop) | peel_sap_stop == "",
      "missing",
      peel_sap_stop
    )
  ) %>%
  count(
    outlier_method,
    peel_sap_stop,
    name = "n"
  ) %>%
  group_by(outlier_method) %>%
  mutate(
    proportion = n / sum(n)
  ) %>%
  ungroup() %>%
  add_outlier_display()


stop_levels <- sort(unique(sap_stop_data$peel_sap_stop))

stop_color_values <- c(
  COL_BLUE_DARK,
  COL_BLUE,
  COL_BLUE_LIGHT,
  COL_ORANGE,
  COL_ORANGE_DARK,
  COL_PURPLE,
  COL_GREY_DARK,
  COL_GREY,
  COL_GREY_LIGHT
)

stop_palette <- setNames(
  rep(stop_color_values, length.out = length(stop_levels)),
  stop_levels
)


figA7_stops <- ggplot(
  sap_stop_data,
  aes(
    x = outlier_label_plot,
    y = proportion,
    fill = peel_sap_stop
  )
) +
  geom_col(
    width = 0.72,
    colour = "white",
    linewidth = 0.2
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    values = stop_palette,
    breaks = stop_levels
  ) +
  labs(
    x = NULL,
    y = "Proportion of Monte Carlo draws",
    fill = "Stopping reason"
  ) +
  theme_paper() +
  theme(
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )


save_plot(
  figA7_stops,
  file.path(
    fig_supp_dir,
    "04_figA7_sap_stop_reasons.pdf"
  ),
  width = 8.5,
  height = 5.0
)

utils::write.csv(
  sap_stop_data,
  file.path(data_dir, "04_figA7_sap_stop_reason_data.csv"),
  row.names = FALSE
)


# ==============================================================================
# 18. Supplementary Figure A8: SAP bias advantage heatmap
# ==============================================================================

sap_advantage_heatmap <- estimation_cell %>%
  filter(estimator_id %in% c("full", "mis_sap")) %>%
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
  mutate(
    sap_bias_advantage = full - mis_sap,
    x_label = factor(
      unname(x_labels_plot[as.character(x_type)]),
      levels = unname(x_labels_plot[x_order])
    ),
    error_label = factor(
      unname(error_labels_plot[as.character(error_type)]),
      levels = unname(error_labels_plot[error_order])
    ),
    cell_label = ifelse(
      is.finite(sap_bias_advantage),
      sprintf("%+.3f", sap_bias_advantage),
      ""
    )
  ) %>%
  add_outlier_display()


max_abs_advantage <- max(
  abs(sap_advantage_heatmap$sap_bias_advantage),
  na.rm = TRUE
)

if (!is.finite(max_abs_advantage) || max_abs_advantage <= 0) {
  max_abs_advantage <- 1
}


figA8_advantage <- ggplot(
  sap_advantage_heatmap,
  aes(
    x = x_label,
    y = error_label,
    fill = sap_bias_advantage
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
    ncol = 2
  ) +
  scale_fill_gradient2(
    low = COL_ORANGE,
    mid = COL_NEUTRAL,
    high = COL_BLUE,
    midpoint = 0,
    limits = c(-max_abs_advantage, max_abs_advantage),
    name = "OLS MAB -\nSAP MAB"
  ) +
  labs(
    x = "Predictor distribution",
    y = "Error distribution"
  ) +
  theme_heatmap() +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    )
  )


save_plot(
  figA8_advantage,
  file.path(
    fig_supp_dir,
    "04_figA8_sap_bias_advantage_heatmap.pdf"
  ),
  width = 9.4,
  height = 7.0
)

utils::write.csv(
  sap_advantage_heatmap,
  file.path(data_dir, "04_figA8_sap_bias_advantage_data.csv"),
  row.names = FALSE
)


# ==============================================================================
# 19. Supplementary tables
# ==============================================================================

# A1: mean runtime.
tabA1_runtime <- make_metric_wide_table(
  data = runtime_broad,
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


# A3: SAP by error distribution.
sap_by_error <- sap_cell %>%
  group_by(error_type, outlier_method) %>%
  summarise(
    detection_rate = safe_mean(detection_rate),
    mean_selected_k = safe_mean(mean_selected_k),
    exact_k_rate = safe_mean(exact_k_rate),
    mean_overlap = safe_mean(mean_overlap),
    mean_abs_bias_sap = safe_mean(mean_abs_bias_sap),
    coverage_sap = safe_mean(coverage_sap),
    error_rate = safe_mean(error_rate),
    .groups = "drop"
  ) %>%
  filter(outlier_method != "none") %>%
  mutate(
    Error = factor(
      unname(error_labels_table[as.character(error_type)]),
      levels = unname(error_labels_table[error_order])
    ),
    Scenario = factor(
      unname(outlier_labels_table[as.character(outlier_method)]),
      levels = unname(outlier_labels_table[outlier_order[-1L]])
    )
  ) %>%
  arrange(Scenario, Error)


tabA3_sap_error <- sap_by_error %>%
  transmute(
    Scenario = as.character(Scenario),
    Error = as.character(Error),
    Detection = fmt_pct(detection_rate, digits = 1L),
    `Mean k` = fmt_num(mean_selected_k, digits = 2L),
    `Exact k` = fmt_pct(exact_k_rate, digits = 1L),
    Overlap = fmt_pct(mean_overlap, digits = 1L),
    `Mean abs. bias` = fmt_num(mean_abs_bias_sap, digits = 3L),
    Coverage = fmt_pct(coverage_sap, digits = 1L),
    `Error rate` = fmt_pct(error_rate, digits = 1L)
  )


write_tex_table(
  data = tabA3_sap_error,
  tex_path = file.path(
    tab_supp_dir,
    "04_tabA3_sap_by_error_distribution.tex"
  ),
  caption = paste0(
    "MIS-SAP performance by error distribution and contamination mechanism. ",
    "Entries are arithmetic means across predictor-distribution cells."
  ),
  label = "tab:robust-sap-error-distribution",
  resize_width = TABLE_WIDTH_WIDE,
  align = "llrrrrrrr",
  font_command = "\\scriptsize"
)


# A4: SAP stopping reasons.
tabA4_stops <- sap_stop_data %>%
  mutate(
    Scenario = as.character(outlier_label_table),
    `Stopping reason` = peel_sap_stop,
    Proportion = fmt_pct(proportion, digits = 1L)
  ) %>%
  select(
    Scenario,
    `Stopping reason`,
    n,
    Proportion
  ) %>%
  rename(Count = n)


write_tex_table(
  data = tabA4_stops,
  tex_path = file.path(
    tab_supp_dir,
    "04_tabA4_sap_stop_reasons.tex"
  ),
  caption = "Distribution of MIS-SAP stopping reasons by contamination mechanism.",
  label = "tab:robust-sap-stop-reasons",
  resize_width = TABLE_WIDTH_COMPACT,
  align = "llrr"
)


# ==============================================================================
# 20. Additional machine-readable full results
# ==============================================================================

# Complete estimator-cell table.
complete_cell_table <- estimation_cell %>%
  arrange(
    outlier_method,
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
    "x_type", "error_type", "outlier_method",
    "estimator_id", "estimator_label",
    "estimator_order", "selection_order"
  )
) %>%
  arrange(
    outlier_method,
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
  selection_cell = selection_cell,
  selection_broad = selection_broad,
  overlap_cell = overlap_cell,
  overlap_broad = overlap_broad,
  runtime_cell = runtime_cell,
  runtime_broad = runtime_broad,
  sap_cell = sap_cell,
  sap_broad = sap_broad,
  method_health = method_health,
  sap_stop_data = sap_stop_data,
  figure1_summary = coef_distribution_summary,
  figure2_data = performance_plot_data,
  figure3_summary = selection_distribution_summary,
  figure4_summary = overlap_distribution_summary
)

saveRDS(
  publication_summaries,
  file.path(data_dir, "04_publication_summaries.rds")
)


# ==============================================================================
# 21. Diagnostics, clipping audit, and manifest
# ==============================================================================

input_audit <- data.frame(
  item = c(
    "Rows in primary RDS",
    "Unique x distributions",
    "Unique error distributions",
    "Unique contamination mechanisms",
    "Unique Monte Carlo iteration IDs",
    "Maximum recorded-versus-recomputed bias difference",
    "Coefficient observations outside main-figure x limits",
    "Selected-k observations outside main-figure x limits",
    "Optional summary RDS present",
    "Optional bias-summary RDS present"
  ),
  value = c(
    nrow(sim),
    length(unique(sim$x_type)),
    length(unique(sim$error_type)),
    length(unique(sim$outlier_method)),
    length(unique(sim$iter)),
    max_bias_difference,
    n_coef_outside_main,
    n_k_outside_main,
    file.exists(input_summary_optional),
    file.exists(input_bias_optional)
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  input_audit,
  file.path(diag_dir, "04_input_and_plot_audit.csv"),
  row.names = FALSE
)


clipping_audit <- data.frame(
  figure = c(
    "04_fig1_coefficient_distributions",
    "04_fig3_selected_k_distributions"
  ),
  lower_limit = c(
    coef_limits[[1L]],
    0
  ),
  upper_limit = c(
    coef_limits[[2L]],
    k_upper
  ),
  observations_outside_display = c(
    n_coef_outside_main,
    n_k_outside_main
  ),
  note = c(
    paste0(
      "Main figure uses central ",
      100 * diff(COEF_MAIN_QUANTILES),
      "% display range; full-range supplementary figure is also saved."
    ),
    paste0(
      "Main figure upper limit uses the ",
      100 * K_MAIN_UPPER_QUANTILE,
      "th percentile and includes the true-k reference."
    )
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  clipping_audit,
  file.path(diag_dir, "04_clipping_audit.csv"),
  row.names = FALSE
)


writeLines(
  capture.output(sessionInfo()),
  file.path(diag_dir, "04_session_info.txt")
)


manifest_paths <- list.files(
  output_root,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

manifest <- data.frame(
  relative_path = substring(
    manifest_paths,
    nchar(output_root) + 2L
  ),
  size_bytes = file.info(manifest_paths)$size,
  modified_time = as.character(file.info(manifest_paths)$mtime),
  stringsAsFactors = FALSE
) %>%
  arrange(relative_path)

utils::write.csv(
  manifest,
  file.path(diag_dir, "04_output_manifest.csv"),
  row.names = FALSE
)


cat("\nScript 84 completed successfully.\n")
cat("Primary input:\n  ", input_main, "\n", sep = "")
cat("Output root:\n  ", output_root, "\n", sep = "")
cat(sprintf(
  "Main coefficient figure display limits: [%.4f, %.4f]\n",
  coef_limits[[1L]],
  coef_limits[[2L]]
))
cat(sprintf(
  "Coefficient observations outside main display: %d\n",
  n_coef_outside_main
))
cat(sprintf(
  "Selected-k observations outside main display: %d\n",
  n_k_outside_main
))
cat(
  "All headline summaries use arithmetic means; ",
  "no median-based headline output was generated.\n",
  sep = ""
)
