# ==============================================================================
# File: scripts/85_output_misspecification.R
# Purpose:
#   Generate publication-ready figures and tables for the Script 05 system:
#
#     05   Frozen estimation / robustness experiment
#     05a  Formal independently calibrated misspecification detector experiment
#     05b  Post-processing: size, power, paired exclusivity, boundaries, k effects
#     05c  Correct-model null calibration / specificity check
#
# Publication interpretation:
#   1. Validity: size, calibration, null scaling
#   2. Detection breadth: reach and power-vs-severity
#   3. Generalist vs specialist: exclusive rejection probabilities
#   4. Detection boundary: scenario x n
#   5. Deletion budget k: operating-point selection
#   6. Specificity and invariance / blind set
#   7. Localisation
#   8. Deletion-arm ladder: detector targeting is not a cleaning rule
#   9. Environment heterogeneity
#  10. Frozen-05 estimation appendix
#
# Output structure:
#   output/05_misspecification/
#     figures/main/          Main-paper PDF figures only
#     figures/supplement/    Supplementary PDF figures only
#     tables/main/           Complete LaTeX table[htbp] environments + CSV
#     tables/supplement/     Supplementary table[htbp] environments + CSV
#     data/                  Numeric data behind every exhibit
#     diagnostics/           Session info, checks, notes, manifest
#
# LaTeX requirements:
#   \usepackage{graphicx}
#   \usepackage{booktabs}
#
# HARD CONSTRAINTS:
#   - Publication primary deletion fraction is k/n = 0.025.
#   - 05a and frozen 05 NEVER share a severity axis:
#       05a uses the formal SD-based calibration;
#       frozen 05 used the older MAD-based calibration.
#   - Frozen-05 exhibits never interpret differences across x_type as severity
#     differences. When frozen-05 results are shown, x_type is kept separate.
#   - Frozen-05 material is used only for estimation / deletion-arm quantities.
#   - 05a/05b/05c provide the formal detector evidence.
#   - Endogeneity with no 80% crossing is displayed as a substantive result,
#     not silently dropped.
#   - Generated tables use [htbp]. Figure floats/captions live in the paper,
#     next to the text they support; this script writes figure PDFs only.
#
# Visual design:
#   - No title or subtitle inside plots; captions are handled in LaTeX.
#   - Soft, low-saturation colours.
#   - No red-green contrast.
#   - MIS is consistently soft blue.
# ==============================================================================


# ==============================================================================
# 0. Packages, paths, configuration
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})


resolve_project_root <- function() {
  candidates <- unique(c(
    normalizePath(".", winslash = "/", mustWork = FALSE),
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

# Formal detector inputs: 05a + renewed 05b + 05c
input_05a_summary <- file.path(
  project_root, "output", "05a_detector_summary.rds"
)
input_05a_boundary <- file.path(
  project_root, "output", "05a_detection_boundary.rds"
)
input_05a_null <- file.path(
  project_root, "output", "05a_null_cutoffs.rds"
)
input_05c_specificity <- file.path(
  project_root, "output", "05a_specificity_corrected.rds"
)
input_05c_correct_cutoffs <- file.path(
  project_root, "output", "05a_correct_cutoffs.rds"
)

# Frozen 05: estimation / deletion-arm evidence only
input_05_frozen_summary <- file.path(
  project_root, "output", "05_misspecification_summary.rds"
)

output_root <- file.path(project_root, "output", "05_misspecification")

fig_main_dir <- file.path(output_root, "figures", "main")
fig_supp_dir <- file.path(output_root, "figures", "supplement")
tab_main_dir <- file.path(output_root, "tables", "main")
tab_supp_dir <- file.path(output_root, "tables", "supplement")
data_dir     <- file.path(output_root, "data")
diag_dir     <- file.path(output_root, "diagnostics")

for (path in c(
  fig_main_dir, fig_supp_dir,
  tab_main_dir, tab_supp_dir,
  data_dir, diag_dir
)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

HAS_FROZEN_05 <- file.exists(input_05_frozen_summary)

if (!HAS_FROZEN_05) {
  warning(
    "Frozen Script 05 summary not found: ", input_05_frozen_summary,
    "\nSections 10 and 12 will be skipped."
  )
}


# Publication operating point.
PRIMARY_K <- 0.025

# PDF is the publication format.
SAVE_PNG_PREVIEWS <- FALSE
PNG_DPI <- 320

TABLE_WIDTH_COMPACT <- "0.70\\columnwidth"
TABLE_WIDTH_MEDIUM  <- "0.85\\columnwidth"
TABLE_WIDTH_WIDE    <- "\\columnwidth"

# ------------------------------------------------------------------------------
# Exhibit registry
#
# Reclassifying an exhibit between the body and supplement is one edit here.
# Plot/table construction below refers to registry IDs rather than hard-coded
# main/supplement directories.
# ------------------------------------------------------------------------------

exhibit_registry <- data.frame(
  id = c(
    "fig01_size",
    "fig02_null_scaling",
    "fig03_reach",
    "fig04_power",
    "fig05_generalist_specialist",
    "figS01_exclusive_classical",
    "fig06_boundary",
    "fig07_budget",
    "fig08_specificity",
    "fig09_targeting_localisation",
    "figS02_environment",
    "tab01_size",
    "tabS01_size_error",
    "tabS02_size_x_n",
    "tab02_null_slopes",
    "tab03_reach",
    "tab04_exclusive_specialists",
    "tab_boundary_summary",
    "tab05_budget",
    "tab_budget_by_scenario",
    "tab06_specificity",
    "tab07_invariance",
    "text07_invariance_argument",
    "tab08_localisation",
    "tabS05_power_error",
    "tabS06_power_x",
    "tabS07_power_n",
    "tabS08_worst_environment",
    "tabA01_mm_correct_wrong",
    "tabA02_ols_mm",
    "tabA03_mm_wins",
    "tabS09_correct_cutoff"
  ),
  type = c(
    rep("figure", 11L),
    rep("table", 11L),
    "text",
    rep("table", 9L)
  ),
  filename = c(
    "05_figS01_size_calibration.pdf",
    "05_figS02_null_scaling.pdf",
    "05_fig01_detection_reach.pdf",
    "05_figS03_power_severity.pdf",
    "05_fig02_generalist_specialist.pdf",
    "05_figS04_exclusive_classical.pdf",
    "05_fig03_boundary_by_n.pdf",
    "05_fig04_budget_reach.pdf",
    "05_fig05_specificity.pdf",
    "05_fig06_targeting_localisation.pdf",
    "05_figS05_environment_heterogeneity.pdf",
    "05_tab01_size_by_k.tex",
    "05_tabS01_size_by_error.tex",
    "05_tabS02_size_by_x_n.tex",
    "05_tab02_null_scaling_slopes.tex",
    "05_tab03_detection_reach.tex",
    "05_tab04_exclusive_specialists.tex",
    "05_tabS03_boundary_by_method.tex",
    "05_tab05_budget_summary.tex",
    "05_tabS04_budget_by_scenario.tex",
    "05_tab06_specificity.tex",
    "05_tab07_invariance_blind_set.tex",
    "05_invariance_argument.tex",
    "05_tab08_localisation.tex",
    "05_tabS05_power_by_error.tex",
    "05_tabS06_power_by_x.tex",
    "05_tabS07_power_by_n.tex",
    "05_tabS08_worst_case_environment.tex",
    "05_tabA01_mm_correct_vs_wrong.tex",
    "05_tabA02_ols_vs_mm.tex",
    "05_tabA03_mm_win_counts.tex",
    "05_tabS09_correct_cutoff_audit.tex"
  ),
  destination = c(
    "supplement",
    "supplement",
    "main",
    "supplement",
    "main",
    "supplement",
    "main",
    "main",
    "main",
    "main",
    "supplement",
    "main",
    "supplement",
    "supplement",
    "main",
    "main",
    "main",
    "supplement",
    "main",
    "supplement",
    "main",
    "main",
    "main",
    "main",
    "supplement",
    "supplement",
    "supplement",
    "supplement",
    "supplement",
    "supplement",
    "supplement",
    "supplement"
  ),
  description = c(
    "Empirical size and calibration visual; supplementary to the body size table.",
    "Null MIS scaling with sample size; supplementary visual for the body slope table.",
    "Number of scenario-environment combinations attaining 80% power.",
    "Six-method power curves over scenario-specific severity; supplementary operating-characteristic detail.",
    "MIS-only versus specialist-only rejection probabilities.",
    "MIS-only versus classical-only rejection probabilities.",
    "MIS 80% boundary by scenario and sample size.",
    "Deletion-budget reach heatmap.",
    "Wrong-model power versus correct-model false-positive rate.",
    "Localisation paired with the frozen-05 deletion-arm ladder.",
    "MIS power heterogeneity over X, error distribution, and n.",
    "Empirical size by deletion fraction.",
    "Empirical size by error distribution.",
    "MIS size by X distribution and n.",
    "Null log-log scaling slopes.",
    "Eighty-percent detection reach.",
    "MIS versus RESET/BP exclusive rejection.",
    "Boundary reach and conditional median by method.",
    "Detection reach by deletion fraction.",
    "MIS reach and boundary by scenario and k.",
    "Wrong-versus-correct specification rejection.",
    "Invariance and blind-set summary.",
    "Algebraic fitted-span invariance argument.",
    "Localisation precision/lift/recall efficiency.",
    "Power by error distribution.",
    "Power by X distribution.",
    "Power by sample size.",
    "Across-environment MIS power distribution at maximum simulated severity.",
    "Frozen-05 correct-versus-wrong full-MM estimation.",
    "Frozen-05 full OLS versus full MM.",
    "Frozen-05 MM win counts.",
    "Correct-model cutoff audit."
  ),
  stringsAsFactors = FALSE
)

# Revision: body tables carry censoring-aware medians; appendix tables expose
# reacher means/medians so the conditional-on-detection comparison is auditable.
exhibit_registry <- bind_rows(
  exhibit_registry,
  data.frame(
    id = c(
      "tabA04_boundary_mean",
      "tabA05_budget_mean",
      "tabA06_frozen_rmse_mean"
    ),
    type = rep("table", 3L),
    filename = c(
      "05_tabA04_boundary_mean.tex",
      "05_tabA05_budget_mean.tex",
      "05_tabA06_frozen_rmse_mean.tex"
    ),
    destination = rep("supplement", 3L),
    description = c(
      "Boundary censoring audit: reach, censored median, reacher mean and reacher median.",
      "Deletion-budget censoring audit with reacher mean/median alongside the censored median.",
      "Frozen-05 cell-level RMSE mean and median side by side."
    ),
    stringsAsFactors = FALSE
  )
)

# Supplementary version of the boundary figure: retain the original 2x4 facets.
exhibit_registry <- bind_rows(
  exhibit_registry,
  data.frame(
    id = "figS_boundary_facet",
    type = "figure",
    filename = "05_figS06_boundary_by_n_faceted.pdf",
    destination = "supplement",
    description = paste0(
      "Faceted 2x4 version of the MIS 80-percent detection boundary by ",
      "scenario and sample size."
    ),
    stringsAsFactors = FALSE
  )
)

# The censoring-aware scenario summaries are body tables. Their mean-based
# counterparts are the A4/A5 registry entries above.
exhibit_registry$destination[
  exhibit_registry$id %in% c("tab_boundary_summary", "tab_budget_by_scenario")
] <- "main"
exhibit_registry$filename[exhibit_registry$id == "tab_boundary_summary"] <-
  "05_tab_boundary_summary.tex"
exhibit_registry$description[exhibit_registry$id == "tab_boundary_summary"] <-
  "Scenario-level MIS boundary summary using the censoring-aware median."
exhibit_registry$filename[exhibit_registry$id == "tab_budget_by_scenario"] <-
  "05_tab_budget_by_scenario.tex"
exhibit_registry$description[exhibit_registry$id == "tab_budget_by_scenario"] <-
  "Scenario-by-k MIS budget summary using the censoring-aware median."

# Figure descriptions after the presentation revision.
exhibit_registry$filename[exhibit_registry$id == "fig01_size"] <-
  "05_figS01_empirical_size_histogram.pdf"
exhibit_registry$description[exhibit_registry$id == "fig01_size"] <-
  "Distribution of cell-level empirical MIS size under the independently calibrated null."
exhibit_registry$description[exhibit_registry$id == "fig08_specificity"] <-
  "Wrong-model rejection versus correct-model false-positive rates."

if (anyDuplicated(exhibit_registry$id)) {
  stop("Exhibit registry contains duplicate IDs.")
}

exhibit_path <- function(id) {
  row <- exhibit_registry[exhibit_registry$id == id, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("Unknown or duplicated exhibit id: ", id)
  } 
  root <- if (row$type == "figure") {
    if (row$destination == "main") fig_main_dir else fig_supp_dir
  } else {
    if (row$destination == "main") tab_main_dir else tab_supp_dir
  }
  
  file.path(root, row$filename)
}

normalize_soft <- function(x) {
  normalizePath(x, winslash = "/", mustWork = FALSE)
}

all_exhibit_paths <- normalize_soft(
  vapply(
    exhibit_registry$id,
    exhibit_path,
    character(1)
  )
)

if (anyDuplicated(all_exhibit_paths)) {
  stop("Two exhibit IDs resolve to one path.")
}

# ------------------------------------------------------------------------------
# Soft colour palette
# ------------------------------------------------------------------------------

COL_MIS       <- "#7FA6C9"  # soft blue
COL_COOK      <- "#D8AE86"  # muted peach
COL_DFBETAS   <- "#B3A7D1"  # soft lavender
COL_LEVERAGE  <- "#B6BCC2"  # cool grey
COL_RESET     <- "#D7A0AD"  # dusty rose
COL_BP        <- "#93B9A4"  # muted sage

COL_BLUE_DARK <- "#527A9D"
COL_BLUE_PALE <- "#DCE8F2"
COL_ROSE_PALE <- "#F0DCE1"
COL_SAGE_PALE <- "#DDE9E0"
COL_PEACH     <- "#E9C9AA"
COL_LILAC     <- "#DDD6EA"
COL_GREY      <- "#7C848B"
COL_GREY_DARK <- "#4E555B"
COL_GREY_PALE <- "#ECEFF1"
COL_NEUTRAL   <- "#F7F5F2"
COL_BLACK     <- "#1F2529"

METHOD_COLOURS <- c(
  "MIS" = COL_MIS,
  "Cook's D" = COL_COOK,
  "DFBETAS" = COL_DFBETAS,
  "Leverage" = COL_LEVERAGE,
  "RESET" = COL_RESET,
  "BP" = COL_BP
)

METHOD_LINETYPES <- c(
  "MIS" = "solid",
  "Cook's D" = "22",
  "DFBETAS" = "42",
  "Leverage" = "13",
  "RESET" = "longdash",
  "BP" = "twodash"
)


# ==============================================================================
# 1. Shared helper functions
# ==============================================================================

near_num <- function(x, target, tol = 1e-12) {
  is.finite(x) & abs(x - target) < tol
}


safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}


safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  stats::median(x)
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


# Boundary summaries must respect right-censoring. A missing boundary means the
# 80% crossing lies beyond the largest simulated severity, not that the cell
# should disappear from the distribution. The median is identified whenever
# fewer than half of environments are censored. Reacher-only means/medians are
# retained for appendix sensitivity tables, but are not the primary summaries.
boundary_stats <- function(boundary, censor_limit) {
  if (length(censor_limit) != 1L || !is.finite(censor_limit)) {
    stop("boundary_stats() requires one finite censor_limit.")
  }
  
  total <- length(boundary)
  reacher_values <- boundary[is.finite(boundary)]
  reachers <- length(reacher_values)
  censored <- total - reachers
  median_identified <- total > 0L && censored < total / 2
  
  censored_median <- if (median_identified) {
    stats::median(c(
      reacher_values,
      rep(censor_limit, censored)
    ))
  } else {
    NA_real_
  }
  
  data.frame(
    reachers = reachers,
    total = total,
    censored = censored,
    reach_rate = if (total > 0L) reachers / total else NA_real_,
    censor_limit = censor_limit,
    median_identified = median_identified,
    censored_median = censored_median,
    reacher_mean = safe_mean(reacher_values),
    reacher_median = safe_median(reacher_values),
    stringsAsFactors = FALSE
  )
}


fmt_num <- function(x, digits = 3L) {
  ifelse(
    is.finite(x),
    formatC(x, format = "f", digits = digits),
    "--"
  )
}


fmt_censored_boundary <- function(
    censored_median,
    median_identified,
    censor_limit,
    digits = 3L
) {
  ifelse(
    median_identified & is.finite(censored_median),
    formatC(censored_median, format = "f", digits = digits),
    paste0(">", formatC(censor_limit, format = "f", digits = digits))
  )
}


fmt_pct <- function(x, digits = 1L) {
  ifelse(
    is.finite(x),
    paste0(formatC(100 * x, format = "f", digits = digits), "%"),
    "--"
  )
}


fmt_sci <- function(x, digits = 2L) {
  ifelse(
    is.finite(x),
    formatC(x, format = "e", digits = digits),
    "--"
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

.label_registry <- character(0)

write_tex_table <- function(
    data,
    tex_path,
    caption,
    label,
    resize_width = TABLE_WIDTH_COMPACT,
    align = NULL,
    csv_path = sub("\\.tex$", ".csv", tex_path),
    placement = "htbp",
    note = NULL,
    note_raw = FALSE,
    addlinespace_after = integer(0)
) {
  if (!is.data.frame(data) || ncol(data) == 0L || nrow(data) == 0L) {
    stop("write_tex_table() requires a non-empty data.frame.")
  }
  
  if (label %in% .label_registry) {
    stop("Duplicate LaTeX label: ", label)
  }
  
  .label_registry <<- c(.label_registry, label)
  
  utils::write.csv(data, csv_path, row.names = FALSE, na = "")
  
  if (is.null(align)) {
    align <- paste0(
      "l",
      paste(rep("r", max(0L, ncol(data) - 1L)), collapse = "")
    )
  }
  
  if (nchar(align) != ncol(data)) {
    stop("Length of LaTeX alignment string must equal number of columns.")
  }
  
  escaped_names <- escape_latex(names(data))
  escaped_data <- lapply(data, escape_latex)
  escaped_data <- as.data.frame(escaped_data, stringsAsFactors = FALSE)
  
  body_rows <- apply(escaped_data, 1L, function(row) {
    paste0(paste(row, collapse = " & "), " \\\\")
  })
  
  if (length(addlinespace_after) > 0L) {
    body_rows <- unlist(
      lapply(seq_along(body_rows), function(i) {
        c(
          body_rows[[i]],
          if (i %in% addlinespace_after) "\\addlinespace" else NULL
        )
      }),
      use.names = FALSE
    )
  }
  
  note_lines <- character(0)
  if (!is.null(note) && length(note) == 1L && nzchar(note)) {
    note_text <- if (isTRUE(note_raw)) note else escape_latex(note)
    note_lines <- c(
      "\\par\\smallskip",
      "\\begin{minipage}{0.96\\linewidth}",
      "\\footnotesize",
      paste0("\\textit{Note:} ", note_text),
      "\\end{minipage}"
    )
  }
  
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
    note_lines,
    "\\end{table}"
  )
  
  writeLines(latex, tex_path, useBytes = TRUE)
  invisible(tex_path)
}


save_plot <- function(plot, pdf_path, width, height) {
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
    png_path <- sub("\\.pdf$", ".png", pdf_path)
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


save_stacked_plots <- function(
    top_plot,
    bottom_plot,
    pdf_path,
    width,
    height,
    top_share = 0.38
) {
  
  draw_stacked <- function() {
    grid::grid.newpage()
    
    layout <- grid::grid.layout(
      nrow = 2L,
      ncol = 1L,
      heights = grid::unit(
        c(top_share, 1 - top_share),
        "null"
      )
    )
    
    grid::pushViewport(
      grid::viewport(layout = layout)
    )
    
    print(
      top_plot,
      vp = grid::viewport(
        layout.pos.row = 1L,
        layout.pos.col = 1L
      )
    )
    
    print(
      bottom_plot,
      vp = grid::viewport(
        layout.pos.row = 2L,
        layout.pos.col = 1L
      )
    )
    
    grid::popViewport()
  }
  
  
  grDevices::pdf(
    pdf_path,
    width = width,
    height = height,
    onefile = FALSE
  )
  
  tryCatch(
    draw_stacked(),
    finally = grDevices::dev.off()
  )
  
  
  if (isTRUE(SAVE_PNG_PREVIEWS)) {
    
    png_path <- sub(
      "\\.pdf$",
      ".png",
      pdf_path
    )
    
    grDevices::png(
      filename = png_path,
      width = width,
      height = height,
      units = "in",
      res = PNG_DPI
    )
    
    tryCatch(
      draw_stacked(),
      finally = grDevices::dev.off()
    )
  }
  
  invisible(pdf_path)
}


theme_85 <- function(base_size = 10.5, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = base_size + 0.5, colour = COL_BLACK),
      axis.text = element_text(size = base_size, colour = COL_GREY_DARK),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.25),
      legend.position = "bottom",
      legend.box = "horizontal",
      strip.background = element_rect(
        fill = COL_NEUTRAL,
        colour = NA
      ),
      strip.text = element_text(
        size = base_size,
        face = "bold",
        colour = COL_GREY_DARK
      ),
      panel.spacing = grid::unit(1.1, "lines"),
      plot.margin = margin(8, 12, 10, 8)
    )
}


theme_85_heatmap <- function(base_size = 10.5, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = base_size + 0.5, colour = COL_BLACK),
      axis.text = element_text(size = base_size, colour = COL_GREY_DARK),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.25),
      legend.position = "bottom",
      strip.background = element_rect(fill = COL_NEUTRAL, colour = NA),
      strip.text = element_text(size = base_size, face = "bold"),
      panel.grid = element_blank(),
      panel.spacing = grid::unit(1.0, "lines"),
      plot.margin = margin(8, 12, 10, 8)
    )
}


# ==============================================================================
# 2. Labels and load/validate inputs
# ==============================================================================

scenario_levels <- c(
  "endogeneity",
  "endogeneity_nl",
  "heterogeneous",
  "heteroskedastic",
  "missing_interaction",
  "nonlinear",
  "ovb",
  "threshold"
)

scenario_labels <- c(
  "endogeneity" = "Linear endogeneity",
  "endogeneity_nl" = "Nonlinear confounding",
  "heterogeneous" = "Heterogeneous slope",
  "heteroskedastic" = "Heteroskedasticity",
  "missing_interaction" = "Missing interaction",
  "nonlinear" = "Nonlinearity",
  "ovb" = "Omitted variable",
  "threshold" = "Threshold"
)

frozen_scenario_labels <- c(
  scenario_labels,
  "structural_break" = "Structural break"
)

x_labels <- c(
  "normal" = "Normal X",
  "mixed_normal" = "Mixed-normal X",
  "contaminated" = "Contaminated X"
)

error_labels <- c(
  "normal" = "Normal",
  "mixed_normal" = "Mixed normal",
  "skewed_t" = "Skewed t",
  "golm" = "GOLM",
  "beta_logistic" = "Beta-logistic",
  "gpd" = "GPD",
  "contaminated" = "Contaminated",
  "pareto" = "Pareto"
)


S <- readRDS(input_05a_summary)
BD <- readRDS(input_05a_boundary)
N0 <- readRDS(input_05a_null)
SPEC <- readRDS(input_05c_specificity)
CORRECT_CUT <- readRDS(input_05c_correct_cutoffs)

W <- S %>%
  filter(model_state == "wrong")

P <- W %>%
  filter(near_num(k_fraction, PRIMARY_K))

PW <- P %>%
  filter(severity_index > 1L)

BD_PRIMARY <- BD %>%
  filter(near_num(k_fraction, PRIMARY_K))

cat(
  "Loaded formal detector summary: ",
  format(nrow(S), big.mark = ","), " rows\n",
  "Loaded renewed boundary table: ",
  format(nrow(BD), big.mark = ","), " rows\n",
  "Primary publication k/n: ",
  PRIMARY_K, "\n",
  sep = ""
)


# ==============================================================================
# 3. SIZE & CALIBRATION
# ==============================================================================

size_raw <- W %>%
  filter(severity_index == 1L)

size_by_k <- size_raw %>%
  group_by(k_fraction) %>%
  summarise(
    MIS = safe_mean(rej_MIS),
    Cook = safe_mean(rej_Cook),
    DFBETAS = safe_mean(rej_DFBETAS),
    Leverage = safe_mean(rej_Leverage),
    RESET_cal = safe_mean(rej_RESET_cal),
    BP_cal = safe_mean(rej_BP_cal),
    RESET_nom = safe_mean(rej_RESET_nom),
    BP_nom = safe_mean(rej_BP_nom),
    .groups = "drop"
  ) %>%
  arrange(k_fraction)

utils::write.csv(
  size_by_k,
  file.path(data_dir, "05_size_by_k.csv"),
  row.names = FALSE
)

# Distribution of empirical MIS size across the 768 severity-0 cells at the
# publication k. Independent calibration should preserve an average near 5%
# without mechanically forcing every evaluation cell to equal 0.05.
size_primary_hist <- size_raw %>%
  filter(near_num(k_fraction, PRIMARY_K))

size_eval_draws <- if ("draws" %in% names(size_primary_hist)) {
  as.integer(round(stats::median(size_primary_hist$draws, na.rm = TRUE)))
} else {
  500L
}

size_target <- 0.05
size_binomial_se <- sqrt(
  size_target * (1 - size_target) / size_eval_draws
)
size_pooled_mean <- safe_mean(size_primary_hist$rej_MIS)
size_observed_sd <- stats::sd(size_primary_hist$rej_MIS, na.rm = TRUE)
size_extra_sd <- sqrt(max(
  size_observed_sd^2 - size_binomial_se^2,
  0
))

size_hist_summary <- data.frame(
  cells = nrow(size_primary_hist),
  evaluation_draws_per_cell = size_eval_draws,
  pooled_mean = size_pooled_mean,
  observed_sd = size_observed_sd,
  binomial_se_at_005 = size_binomial_se,
  residual_dispersion = size_extra_sd,
  lower_two_mcse = size_target - 2 * size_binomial_se,
  upper_two_mcse = size_target + 2 * size_binomial_se,
  stringsAsFactors = FALSE
)

utils::write.csv(
  size_primary_hist %>%
    select(
      cell_id, env_id, n, x_type, error_type, scenario,
      k_fraction, rej_MIS
    ),
  file.path(data_dir, "05_empirical_size_histogram_primary_k.csv"),
  row.names = FALSE
)

utils::write.csv(
  size_hist_summary,
  file.path(data_dir, "05_empirical_size_histogram_summary.csv"),
  row.names = FALSE
)

fig_size <- ggplot(
  size_primary_hist,
  aes(x = rej_MIS)
) +
  geom_histogram(
    bins = 30,
    fill = COL_MIS,
    colour = "white",
    linewidth = 0.35,
    alpha = 0.92
  ) +
  geom_vline(
    xintercept = size_target,
    colour = COL_GREY_DARK,
    linewidth = 0.75
  ) +
  geom_vline(
    xintercept = c(
      size_target - 2 * size_binomial_se,
      size_target + 2 * size_binomial_se
    ),
    colour = COL_GREY,
    linetype = "dashed",
    linewidth = 0.65
  ) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  labs(
    x = "Empirical rejection rate at severity 0",
    y = "Cells"
  ) +
  theme_85()

save_plot(
  fig_size,
  exhibit_path("fig01_size"),
  width = 7.8,
  height = 4.9
)


size_table_main <- size_by_k %>%
  transmute(
    `k/n` = fmt_pct(k_fraction, 1),
    MIS = fmt_pct(MIS, 2),
    `Cook's D` = fmt_pct(Cook, 2),
    DFBETAS = fmt_pct(DFBETAS, 2),
    Leverage = fmt_pct(Leverage, 2),
    `RESET (cal.)` = fmt_pct(RESET_cal, 2),
    `RESET (nom.)` = fmt_pct(RESET_nom, 2),
    `BP (cal.)` = fmt_pct(BP_cal, 2),
    `BP (nom.)` = fmt_pct(BP_nom, 2)
  )

write_tex_table(
  size_table_main,
  exhibit_path("tab01_size"),
  caption = paste0(
    "Empirical size by deletion fraction. Detector cutoffs are estimated from ",
    "an independent calibration stream."
  ),
  label = "tab:05-size-by-k",
  resize_width = TABLE_WIDTH_WIDE,
  align = "lrrrrrrrr",
  placement = "htbp"
)


size_primary <- size_raw %>%
  filter(near_num(k_fraction, PRIMARY_K))

size_by_error <- size_primary %>%
  group_by(error_type) %>%
  summarise(
    MIS = safe_mean(rej_MIS),
    RESET_cal = safe_mean(rej_RESET_cal),
    RESET_nom = safe_mean(rej_RESET_nom),
    BP_cal = safe_mean(rej_BP_cal),
    BP_nom = safe_mean(rej_BP_nom),
    .groups = "drop"
  )

utils::write.csv(
  size_by_error,
  file.path(data_dir, "05_size_by_error_primary_k.csv"),
  row.names = FALSE
)

write_tex_table(
  size_by_error %>%
    transmute(
      `Error distribution` = unname(error_labels[error_type]),
      MIS = fmt_pct(MIS, 2),
      `RESET (cal.)` = fmt_pct(RESET_cal, 2),
      `RESET (nom.)` = fmt_pct(RESET_nom, 2),
      `BP (cal.)` = fmt_pct(BP_cal, 2),
      `BP (nom.)` = fmt_pct(BP_nom, 2)
    ),
  exhibit_path("tabS01_size_error"),
  caption = paste0(
    "Empirical size by error distribution at k/n = 2.5 percent for MIS, ",
    "RESET, and BP."
  ),
  label = "tab:05-size-by-error",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrrr",
  placement = "htbp"
)


size_by_x_n <- size_primary %>%
  group_by(x_type, n) %>%
  summarise(
    MIS = safe_mean(rej_MIS),
    .groups = "drop"
  )

utils::write.csv(
  size_by_x_n,
  file.path(data_dir, "05_size_by_x_n_primary_k.csv"),
  row.names = FALSE
)

write_tex_table(
  size_by_x_n %>%
    mutate(
      x_order = factor(
        x_type,
        levels = c(
          "normal",
          "mixed_normal",
          "contaminated"
        )
      )
    ) %>%
    arrange(x_order, n) %>%
    transmute(
      `X distribution` = unname(x_labels[x_type]),
      n = n,
      `MIS size` = fmt_pct(MIS, 2)
    ),
  exhibit_path("tabS02_size_x_n"),
  caption = "MIS empirical size by X distribution and sample size at k/n = 2.5 percent.",
  label = "tab:05-size-by-x-n",
  resize_width = TABLE_WIDTH_COMPACT,
  align = "lrr",
  placement = "htbp",
  addlinespace_after = c(4L, 8L)
)


# ------------------------------------------------------------------------------
# Null scaling
# ------------------------------------------------------------------------------

null_scale <- N0 %>%
  filter(near_num(k_fraction, PRIMARY_K)) %>%
  group_by(x_type, n) %>%
  summarise(
    null_mean_MIS = safe_mean(null_mean_MIS),
    .groups = "drop"
  ) %>%
  arrange(x_type, n)

null_slopes <- bind_rows(
  lapply(
    split(null_scale, null_scale$x_type),
    function(g) {
      fit <- stats::lm(log(null_mean_MIS) ~ log(n), data = g)
      data.frame(
        x_type = g$x_type[[1L]],
        slope = unname(stats::coef(fit)[[2L]]),
        r_squared = summary(fit)$r.squared,
        stringsAsFactors = FALSE
      )
    }
  )
)

null_scale <- null_scale %>%
  left_join(null_slopes, by = "x_type") %>%
  mutate(
    x_display = paste0(
      unname(x_labels[x_type]),
      " (slope ", formatC(slope, format = "f", digits = 3), ")"
    )
  )

utils::write.csv(
  null_scale,
  file.path(data_dir, "05_null_scaling.csv"),
  row.names = FALSE
)

utils::write.csv(
  null_slopes,
  file.path(data_dir, "05_null_scaling_slopes.csv"),
  row.names = FALSE
)

fig_null_scale <- ggplot(
  null_scale,
  aes(
    x = n,
    y = null_mean_MIS,
    colour = x_display,
    group = x_display
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.3) +
  scale_x_log10(
    breaks = sort(unique(null_scale$n)),
    labels = scales::comma
  ) +
  scale_y_log10() +
  scale_colour_manual(
    values = c(
      "#7FA6C9",
      "#A9BFCF",
      "#B3A7D1"
    )
  ) +
  labs(
    x = "Sample size n (log scale)",
    y = "Mean null MIS statistic (log scale)",
    colour = NULL
  ) +
  theme_85()

save_plot(
  fig_null_scale,
  exhibit_path("fig02_null_scaling"),
  width = 7.2,
  height = 4.8
)


write_tex_table(
  null_slopes %>%
    mutate(
      x_order = factor(
        x_type,
        levels = c("normal", "mixed_normal", "contaminated")
      )
    ) %>%
    arrange(x_order) %>%
    transmute(
      `X distribution` = unname(x_labels[x_type]),
      `Log-log slope` = fmt_num(slope, 3),
      `R-squared` = fmt_num(r_squared, 4)
    ),
  exhibit_path("tab02_null_slopes"),
  caption = "Log-log scaling of the null mean MIS statistic with sample size.",
  label = "tab:05-null-scaling",
  resize_width = TABLE_WIDTH_COMPACT,
  align = "lrr",
  placement = "htbp"
)


# ==============================================================================
# 4. DETECTION BREADTH
# ==============================================================================

boundary_columns <- c(
  "MIS" = "b_MIS",
  "Cook's D" = "b_Cook",
  "DFBETAS" = "b_DFBETAS",
  "Leverage" = "b_Leverage",
  "RESET" = "b_RESET",
  "BP" = "b_BP"
)

reach_counts <- bind_rows(
  lapply(
    names(boundary_columns),
    function(method) {
      v <- BD_PRIMARY[[boundary_columns[[method]]]]
      data.frame(
        method = method,
        reached = sum(is.finite(v)),
        total = length(v),
        reach_rate = mean(is.finite(v)),
        stringsAsFactors = FALSE
      )
    }
  )
) %>%
  arrange(desc(reached)) %>%
  mutate(
    method = factor(method, levels = rev(method))
  )

utils::write.csv(
  reach_counts,
  file.path(data_dir, "05_reach_counts_primary_k.csv"),
  row.names = FALSE
)

fig_reach <- ggplot(
  reach_counts,
  aes(x = reached, y = method, fill = method)
) +
  geom_col(width = 0.68, alpha = 0.88) +
  geom_text(
    aes(label = paste0(reached, " / ", total)),
    hjust = -0.10,
    size = 3.4,
    colour = COL_GREY_DARK
  ) +
  scale_fill_manual(values = METHOD_COLOURS, guide = "none") +
  scale_x_continuous(
    limits = c(0, max(reach_counts$total) * 1.10),
    breaks = seq(0, max(reach_counts$total), by = 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Environments attaining 80% power",
    y = NULL
  ) +
  theme_85() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

save_plot(
  fig_reach,
  exhibit_path("fig03_reach"),
  width = 7.1,
  height = 4.7
)


write_tex_table(
  reach_counts %>%
    arrange(desc(reached)) %>%
    transmute(
      Method = as.character(method),
      `Reached 80 percent` = reached,
      Total = total,
      `Reach rate` = fmt_pct(reach_rate, 1)
    ),
  exhibit_path("tab03_reach"),
  caption = "Eighty-percent detection reach at k/n = 2.5 percent.",
  label = "tab:05-detection-reach",
  resize_width = TABLE_WIDTH_COMPACT,
  align = "lrrr",
  placement = "htbp"
)


# ------------------------------------------------------------------------------
# Power-vs-severity curves
# ------------------------------------------------------------------------------

power_metric_lookup <- c(
  "rej_MIS" = "MIS",
  "rej_Cook" = "Cook's D",
  "rej_DFBETAS" = "DFBETAS",
  "rej_Leverage" = "Leverage",
  "rej_RESET_cal" = "RESET",
  "rej_BP_cal" = "BP"
)

power_curve <- P %>%
  select(
    scenario,
    severity_index,
    severity_target,
    all_of(names(power_metric_lookup))
  ) %>%
  pivot_longer(
    cols = all_of(names(power_metric_lookup)),
    names_to = "metric",
    values_to = "rejection_rate"
  ) %>%
  mutate(
    method = unname(power_metric_lookup[metric])
  ) %>%
  group_by(
    scenario,
    severity_index,
    severity_target,
    method
  ) %>%
  summarise(
    power = safe_mean(rejection_rate),
    .groups = "drop"
  ) %>%
  mutate(
    scenario_display = factor(
      unname(scenario_labels[scenario]),
      levels = unname(scenario_labels[scenario_levels])
    ),
    method = factor(
      method,
      levels = c("MIS", "Cook's D", "DFBETAS", "Leverage", "RESET", "BP")
    )
  )

utils::write.csv(
  power_curve,
  file.path(data_dir, "05_power_severity_curves_primary_k.csv"),
  row.names = FALSE
)

fig_power <- ggplot(
  power_curve,
  aes(
    x = severity_target,
    y = power,
    colour = method,
    linetype = method,
    group = method
  )
) +
  geom_hline(
    yintercept = 0.80,
    linetype = "dotted",
    colour = COL_GREY,
    linewidth = 0.45
  ) +
  geom_line(linewidth = 0.62) +
  geom_point(size = 1.15, alpha = 0.90) +
  facet_wrap(
    ~ scenario_display,
    ncol = 4,
    scales = "free_x",
    drop = FALSE
  ) +
  scale_colour_manual(values = METHOD_COLOURS) +
  scale_linetype_manual(values = METHOD_LINETYPES) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = "Misspecification severity (scenario-specific scale)",
    y = "Empirical power",
    colour = NULL,
    linetype = NULL
  ) +
  theme_85(base_size = 9.4) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8.5),
    strip.text = element_text(size = 9.2, face = "bold")
  )

save_plot(
  fig_power,
  exhibit_path("fig04_power"),
  width = 8.8,
  height = 5.4
)


# ==============================================================================
# 5. GENERALIST VS SPECIALIST: EXCLUSIVE REJECTION
# ==============================================================================

exclusive_pairs <- list(
  "Cook's D" = c("only_MIS_vs_Cook", "only_Cook_vs_MIS"),
  "DFBETAS" = c("only_MIS_vs_DFBETAS", "only_DFBETAS_vs_MIS"),
  "Leverage" = c("only_MIS_vs_Leverage", "only_Leverage_vs_MIS"),
  "RESET" = c("only_MIS_vs_RESET", "only_RESET_vs_MIS"),
  "BP" = c("only_MIS_vs_BP", "only_BP_vs_MIS")
)


exclusive_summary <- bind_rows(
  lapply(
    names(exclusive_pairs),
    function(rival) {
      cols <- exclusive_pairs[[rival]]
      
      PW %>%
        group_by(scenario) %>%
        summarise(
          mis_only = safe_mean(.data[[cols[[1L]]]]),
          rival_only = safe_mean(.data[[cols[[2L]]]]),
          .groups = "drop"
        ) %>%
        mutate(rival = rival)
    }
  )
) %>%
  mutate(
    scenario_display = factor(
      unname(scenario_labels[scenario]),
      levels = rev(unname(scenario_labels[scenario_levels]))
    )
  )

utils::write.csv(
  exclusive_summary,
  file.path(data_dir, "05_exclusive_rejection_primary_k.csv"),
  row.names = FALSE
)


make_exclusive_plot <- function(data, rivals, exhibit_id, width) {
  d <- bind_rows(
    data %>%
      filter(rival %in% rivals) %>%
      transmute(
        scenario,
        scenario_display,
        rival,
        side = "MIS only",
        signed_rate = mis_only
      ),
    data %>%
      filter(rival %in% rivals) %>%
      transmute(
        scenario,
        scenario_display,
        rival,
        side = "Rival only",
        signed_rate = -rival_only
      )
  )
  
  p <- ggplot(
    d,
    aes(
      x = signed_rate,
      y = scenario_display,
      fill = side
    )
  ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.55,
      colour = COL_GREY_DARK
    ) +
    geom_col(
      width = 0.70,
      alpha = 0.90
    ) +
    facet_wrap(
      ~ rival,
      nrow = 1,
      scales = "free_x"
    ) +
    scale_fill_manual(
      values = c(
        "MIS only" = COL_BLUE_PALE,
        "Rival only" = COL_ROSE_PALE
      )
    ) +
    scale_x_continuous(
      labels = function(x) scales::percent(abs(x), accuracy = 1),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    labs(
      x = "Rival only  <-    |    ->  MIS only",
      y = NULL,
      fill = NULL
    ) +
    theme_85(base_size = 9.6) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 9.6, face = "bold"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  save_plot(
    p,
    exhibit_path(exhibit_id),
    width = width,
    height = 5.8
  )
}


make_exclusive_plot(
  exclusive_summary,
  rivals = c("RESET", "BP"),
  exhibit_id = "fig05_generalist_specialist",
  width = 9.0
)


make_exclusive_plot(
  exclusive_summary,
  rivals = c("Cook's D", "DFBETAS", "Leverage"),
  exhibit_id = "figS01_exclusive_classical",
  width = 10.8
)


exclusive_specialist_table <- exclusive_summary %>%
  filter(rival %in% c("RESET", "BP")) %>%
  mutate(
    scenario_order = match(scenario, scenario_levels),
    rival_order = match(rival, c("BP", "RESET"))
  ) %>%
  arrange(rival_order, scenario_order) %>%
  transmute(
    Scenario = unname(scenario_labels[scenario]),
    Rival = rival,
    `MIS only` = fmt_pct(mis_only, 1),
    `Rival only` = fmt_pct(rival_only, 1)
  )

write_tex_table(
  exclusive_specialist_table,
  exhibit_path("tab04_exclusive_specialists"),
  caption = paste0(
    "Exclusive rejection probabilities for MIS versus RESET and BP at ",
    "k/n = 2.5 percent."
  ),
  label = "tab:05-exclusive-specialists",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "llrr",
  placement = "htbp",
  addlinespace_after = 8L
)


# ==============================================================================
# 6. EIGHTY-PERCENT DETECTION BOUNDARY
# ==============================================================================

# Maximum simulated severity is the right-censoring point for a scenario.
scenario_max_severity <- W %>%
  group_by(scenario) %>%
  summarise(
    max_simulated_severity = max(severity_target, na.rm = TRUE),
    .groups = "drop"
  )

BD_CENSOR <- BD %>%
  left_join(scenario_max_severity, by = "scenario")

BD_PRIMARY_CENSOR <- BD_CENSOR %>%
  filter(near_num(k_fraction, PRIMARY_K))


# ------------------------------------------------------------------------------
# 6.1 Figure: censored median by n
# ------------------------------------------------------------------------------

boundary_by_n <- BD_PRIMARY_CENSOR %>%
  group_by(scenario, n, max_simulated_severity) %>%
  group_modify(
    ~ boundary_stats(
      .x$b_MIS,
      .y$max_simulated_severity[[1L]]
    )
  ) %>%
  ungroup() %>%
  mutate(
    scenario_display = factor(
      unname(scenario_labels[scenario]),
      levels = unname(scenario_labels[scenario_levels])
    ),
    n_display = factor(
      n,
      levels = sort(unique(n))
    ),
    plot_boundary = ifelse(
      median_identified,
      censored_median,
      censor_limit
    )
  )

utils::write.csv(
  boundary_by_n,
  file.path(data_dir, "05_boundary_by_n_primary_k.csv"),
  row.names = FALSE
)
# Six scenarios have at least one identified censoring-aware median.
boundary_line_scenarios <- boundary_by_n %>%
  group_by(scenario) %>%
  summarise(
    any_identified = any(median_identified),
    .groups = "drop"
  ) %>%
  filter(any_identified) %>%
  pull(scenario)


boundary_loglog_slopes <- boundary_by_n %>%
  filter(
    scenario %in% boundary_line_scenarios,
    median_identified,
    is.finite(censored_median),
    censored_median > 0
  ) %>%
  group_by(scenario) %>%
  group_modify(
    ~ {
      fit <- lm(
        log(censored_median) ~ log(n),
        data = .x
      )
      
      data.frame(
        n_points = nrow(.x),
        slope = unname(coef(fit)[["log(n)"]]),
        r_squared = summary(fit)$r.squared
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    scenario_display = unname(scenario_labels[scenario])
  )


utils::write.csv(
  boundary_loglog_slopes,
  file.path(data_dir, "05_boundary_loglog_slopes.csv"),
  row.names = FALSE
)


boundary_scenario_colours <- c(
  "endogeneity_nl" = COL_DFBETAS,
  "heterogeneous" = COL_BP,
  "heteroskedastic" = COL_COOK,
  "missing_interaction" = COL_RESET,
  "nonlinear" = COL_MIS,
  "threshold" = COL_BLUE_DARK
)

boundary_scenario_linetypes <- c(
  "endogeneity_nl" = "solid",
  "heterogeneous" = "22",
  "heteroskedastic" = "42",
  "missing_interaction" = "13",
  "nonlinear" = "longdash",
  "threshold" = "twodash"
)

boundary_scenario_shapes <- c(
  "endogeneity_nl" = 15,
  "heterogeneous" = 16,
  "heteroskedastic" = 17,
  "missing_interaction" = 18,
  "nonlinear" = 8,
  "threshold" = 3
)


boundary_plot_data <- boundary_by_n %>%
  filter(scenario %in% boundary_line_scenarios) %>%
  arrange(scenario, n) %>%
  group_by(scenario) %>%
  mutate(
    identified_boundary = ifelse(
      median_identified,
      censored_median,
      NA_real_
    ),
    line_run = cumsum(!median_identified)
  ) %>%
  ungroup()

fig_boundary <- ggplot(
  boundary_plot_data,
  aes(
    x = n,
    colour = scenario,
    linetype = scenario,
    shape = scenario,
    group = scenario
  )
) +
  geom_line(
    data = boundary_plot_data %>%
      filter(median_identified),
    aes(
      y = censored_median,
      group = interaction(scenario, line_run)
    ),
    linewidth = 0.85
  ) +
  geom_point(
    data = boundary_plot_data %>%
      filter(median_identified),
    aes(y = censored_median),
    size = 2.4
  ) +
  geom_point(
    data = boundary_plot_data %>%
      filter(!median_identified),
    aes(y = censor_limit),
    shape = 24,
    size = 2.8,
    fill = "white",
    stroke = 0.7,
    show.legend = FALSE
  ) +
  scale_x_log10(
    breaks = c(500, 1000, 2500, 5000),
    labels = scales::comma
  ) +
  scale_colour_manual(
    values = boundary_scenario_colours,
    breaks = names(boundary_scenario_linetypes),
    labels = unname(
      scenario_labels[names(boundary_scenario_linetypes)]
    ),
    name = NULL
  ) +
  scale_linetype_manual(
    values = boundary_scenario_linetypes,
    breaks = names(boundary_scenario_linetypes),
    labels = unname(
      scenario_labels[names(boundary_scenario_linetypes)]
    ),
    name = NULL
  ) +
  scale_shape_manual(
    values = boundary_scenario_shapes,
    breaks = names(boundary_scenario_linetypes),
    labels = unname(
      scenario_labels[names(boundary_scenario_linetypes)]
    ),
    name = NULL
  ) +
  labs(
    x = "Sample size n (log scale)",
    y = "Censoring-aware median severity at 80% power",
    colour = NULL,
    linetype = NULL,
    caption = paste0(
      "Linear endogeneity and omitted variable never reached ",
      "80% power at any sample size."
    )
  ) +
  theme_85(base_size = 10.2) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.width = grid::unit(1.25, "cm"),
    legend.key.height = grid::unit(0.45, "cm"),
    plot.caption = element_text(
      size = 8.5,
      colour = COL_GREY_DARK,
      hjust = 0
    )
  )

save_plot(
  fig_boundary,
  exhibit_path("fig06_boundary"),
  width = 7.0,
  height = 5.0
)

fig_boundary_facet <- ggplot(
  boundary_by_n,
  aes(
    x = n_display,
    y = plot_boundary,
    group = 1
  )
) +
  geom_hline(
    data = boundary_by_n %>%
      distinct(
        scenario_display,
        censor_limit
      ),
    aes(yintercept = censor_limit),
    inherit.aes = FALSE,
    linetype = "dashed",
    colour = COL_GREY,
    linewidth = 0.5
  ) +
  geom_line(
    data = boundary_by_n %>%
      filter(median_identified),
    linewidth = 0.75,
    colour = COL_MIS,
    na.rm = TRUE
  ) +
  geom_point(
    data = boundary_by_n %>%
      filter(median_identified),
    size = 2.1,
    colour = COL_BLUE_DARK,
    na.rm = TRUE
  ) +
  geom_point(
    data = boundary_by_n %>%
      filter(!median_identified),
    shape = 24,
    size = 2.6,
    fill = COL_NEUTRAL,
    colour = COL_GREY_DARK,
    stroke = 0.55,
    na.rm = TRUE
  ) +
  geom_text(
    data = boundary_by_n %>%
      filter(!median_identified),
    aes(label = ">max"),
    vjust = -0.8,
    size = 2.5,
    colour = COL_GREY_DARK,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ scenario_display,
    ncol = 4,
    scales = "free_y",
    drop = FALSE
  ) +
  expand_limits(y = 0) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.20))
  ) +
  labs(
    x = "Sample size n",
    y = "Censoring-aware median severity at 80% power"
  ) +
  theme_85(base_size = 9.2) +
  theme(
    strip.text = element_text(
      size = 8.9,
      face = "bold"
    )
  )

save_plot(
  fig_boundary_facet,
  exhibit_path("figS_boundary_facet"),
  width = 8.8,
  height = 5.2
)

boundary_loglog_labels <- setNames(
  paste0(
    unname(
      scenario_labels[
        boundary_loglog_slopes$scenario
      ]
    ),
    " (slope ",
    formatC(
      boundary_loglog_slopes$slope,
      format = "f",
      digits = 2
    ),
    ")"
  ),
  boundary_loglog_slopes$scenario
)


fig_boundary_loglog <- ggplot(
  boundary_plot_data,
  aes(
    x = n,
    colour = scenario,
    linetype = scenario,
    group = scenario
  )
) +
  geom_line(
    data = boundary_plot_data %>%
      filter(median_identified),
    aes(
      y = censored_median,
      group = interaction(scenario, line_run)
    ),
    linewidth = 0.85
  ) +
  geom_point(
    data = boundary_plot_data %>%
      filter(median_identified),
    aes(y = censored_median),
    size = 2.3,
    show.legend = FALSE
  ) +
  geom_point(
    data = boundary_plot_data %>%
      filter(!median_identified),
    aes(y = censor_limit),
    shape = 24,
    size = 2.8,
    fill = "white",
    stroke = 0.7,
    show.legend = FALSE
  ) +
  scale_x_log10(
    breaks = c(500, 1000, 2500, 5000),
    labels = scales::comma
  ) +
  scale_y_log10() +
  scale_colour_manual(
    values = boundary_scenario_colours,
    breaks = names(boundary_scenario_colours),
    labels = boundary_loglog_labels[
      names(boundary_scenario_colours)
    ]
  ) +
  scale_linetype_manual(
    values = boundary_scenario_linetypes,
    breaks = names(boundary_scenario_linetypes),
    labels = boundary_loglog_labels[
      names(boundary_scenario_linetypes)
    ]
  ) +
  labs(
    x = "Sample size n (log scale)",
    y = "Censoring-aware median severity at 80% power (log scale)",
    colour = NULL,
    linetype = NULL,
    caption = paste0(
      "Linear endogeneity and omitted variable never reached ",
      "80% power at any sample size."
    )
  ) +
  theme_85(base_size = 10.2) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.width = grid::unit(2.8, "cm"),
    legend.key.height = grid::unit(0.5, "cm"),
    plot.caption = element_text(
      size = 8.5,
      colour = COL_GREY_DARK,
      hjust = 0
    )
  )

save_plot(
  fig_boundary_loglog,
  file.path(diag_dir, "05_boundary_loglog_check.pdf"),
  width = 7.0,
  height = 5.0
)

# ------------------------------------------------------------------------------
# 6.2 Body table: MIS scenario boundary, censoring retained
# ------------------------------------------------------------------------------

boundary_primary_summary <- BD_PRIMARY_CENSOR %>%
  group_by(scenario, max_simulated_severity) %>%
  group_modify(
    ~ boundary_stats(
      .x$b_MIS,
      .y$max_simulated_severity[[1L]]
    )
  ) %>%
  ungroup() %>%
  arrange(factor(scenario, levels = scenario_levels))

utils::write.csv(
  boundary_primary_summary,
  file.path(data_dir, "05_boundary_censoring_summary_primary_k.csv"),
  row.names = FALSE
)

boundary_note <- paste0(
  "Boundaries are right-censored at the maximum simulated severity, so the ",
  "censoring-aware median is reported; it is identified whenever fewer than ",
  "half of environments are censored. Means over the reaching subset, which ",
  "condition on detection succeeding and therefore understate the boundary, ",
  "are given in Table~\\ref{tab:05-boundary-mean}."
)

write_tex_table(
  boundary_primary_summary %>%
    transmute(
      Scenario = unname(scenario_labels[scenario]),
      `Reachers / 96` = paste0(reachers, " / ", total),
      Censored = censored,
      `Reach rate` = fmt_pct(reach_rate, 1),
      `Censoring-aware median` = fmt_censored_boundary(
        censored_median,
        median_identified,
        censor_limit,
        3
      )
    ),
  exhibit_path("tab_boundary_summary"),
  caption = paste0(
    "MIS 80-percent detection boundary by scenario at k/n = 2.5 percent."
  ),
  label = "tab:05-boundary-summary",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrr",
  placement = "htbp",
  note = boundary_note,
  note_raw = TRUE
)


# Appendix A4: same boundary object, now exposing the reacher-only mean/median.
write_tex_table(
  boundary_primary_summary %>%
    transmute(
      Scenario = unname(scenario_labels[scenario]),
      `Reachers / 96` = paste0(reachers, " / ", total),
      Censored = censored,
      `Reach rate` = fmt_pct(reach_rate, 1),
      `Censoring-aware median` = fmt_censored_boundary(
        censored_median,
        median_identified,
        censor_limit,
        3
      ),
      `Reacher mean` = fmt_num(reacher_mean, 3),
      `Reacher median` = fmt_num(reacher_median, 3)
    ),
  exhibit_path("tabA04_boundary_mean"),
  caption = paste0(
    "Boundary censoring audit at k/n = 2.5 percent. Reacher-only means and ",
    "medians condition on the 80-percent crossing being observed; the censored ",
    "median retains the non-reaching environments."
  ),
  label = "tab:05-boundary-mean",
  resize_width = TABLE_WIDTH_WIDE,
  align = "lrrrrrr",
  placement = "htbp"
)


# ==============================================================================
# 7. DELETION BUDGET k
# ==============================================================================

# Normalized boundaries are used only for the pooled cross-scenario boundary
# summary. Scenario tables remain on each scenario's own severity scale.
BD_BUDGET <- BD_CENSOR %>%
  mutate(
    b_MIS_normalized = ifelse(
      is.finite(b_MIS) &
        is.finite(max_simulated_severity) &
        max_simulated_severity > 0,
      b_MIS / max_simulated_severity,
      NA_real_
    )
  )

budget_mis <- BD_BUDGET %>%
  group_by(scenario, k_fraction, max_simulated_severity) %>%
  group_modify(
    ~ boundary_stats(
      .x$b_MIS,
      .y$max_simulated_severity[[1L]]
    )
  ) %>%
  ungroup()

budget_pooled <- BD_BUDGET %>%
  group_by(k_fraction) %>%
  group_modify(
    ~ boundary_stats(
      .x$b_MIS_normalized,
      1
    )
  ) %>%
  ungroup() %>%
  mutate(scenario = "pooled")

budget_heat <- bind_rows(
  budget_pooled %>%
    select(scenario, k_fraction, reach_rate),
  budget_mis %>%
    select(scenario, k_fraction, reach_rate)
) %>%
  mutate(
    scenario_display = ifelse(
      scenario == "pooled",
      "Pooled",
      unname(scenario_labels[scenario])
    )
  ) %>%
  group_by(scenario_display) %>%
  mutate(
    row_best = max(reach_rate, na.rm = TRUE),
    is_best = is.finite(reach_rate) &
      row_best > 0 &
      abs(reach_rate - row_best) < 1e-12,
    cell_label = paste0(
      fmt_pct(reach_rate, 1),
      ifelse(is_best, " *", "")
    )
  ) %>%
  ungroup()

budget_order <- c(
  "Pooled",
  unname(scenario_labels[scenario_levels])
)

budget_heat <- budget_heat %>%
  mutate(
    scenario_display = factor(
      scenario_display,
      levels = rev(budget_order)
    ),
    k_display = factor(
      fmt_pct(k_fraction, 1),
      levels = fmt_pct(sort(unique(k_fraction)), 1)
    )
  )

utils::write.csv(
  budget_heat,
  file.path(data_dir, "05_budget_mis_reach.csv"),
  row.names = FALSE
)

utils::write.csv(
  budget_mis,
  file.path(data_dir, "05_budget_censoring_by_scenario.csv"),
  row.names = FALSE
)

fig_budget <- ggplot(
  budget_heat,
  aes(
    x = k_display,
    y = scenario_display,
    fill = reach_rate
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.9
  ) +
  geom_text(
    aes(label = cell_label),
    size = 3.25,
    colour = COL_GREY_DARK
  ) +
  scale_fill_gradient(
    low = COL_NEUTRAL,
    high = COL_MIS,
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = scales::percent_format(accuracy = 1),
    guide = guide_colourbar(
      barwidth = grid::unit(7, "cm"),
      barheight = grid::unit(0.45, "cm")
    )
  ) +
  labs(
    x = "Deletion fraction k/n",
    y = NULL,
    fill = "80% reach"
  ) +
  theme_85_heatmap() +
  theme(
    axis.text.x = element_text(size = 9.5),
    axis.text.y = element_text(size = 9.5)
  )

save_plot(
  fig_budget,
  exhibit_path("fig07_budget"),
  width = 7.8,
  height = 5.8
)


# All-method pooled reach remains useful because reach is scale-free.
budget_all_methods <- bind_rows(
  lapply(
    sort(unique(BD_BUDGET$k_fraction)),
    function(kf) {
      g <- BD_BUDGET %>% filter(near_num(k_fraction, kf))
      mis_pool <- boundary_stats(g$b_MIS_normalized, 1)
      
      data.frame(
        k_fraction = kf,
        MIS_reach = mean(is.finite(g$b_MIS)),
        MIS_count = sum(is.finite(g$b_MIS)),
        MIS_censored_median_normalized = mis_pool$censored_median,
        Cook_reach = mean(is.finite(g$b_Cook)),
        Cook_count = sum(is.finite(g$b_Cook)),
        DFBETAS_reach = mean(is.finite(g$b_DFBETAS)),
        DFBETAS_count = sum(is.finite(g$b_DFBETAS)),
        Leverage_reach = mean(is.finite(g$b_Leverage)),
        Leverage_count = sum(is.finite(g$b_Leverage)),
        RESET_reach = mean(is.finite(g$b_RESET)),
        RESET_count = sum(is.finite(g$b_RESET)),
        BP_reach = mean(is.finite(g$b_BP)),
        BP_count = sum(is.finite(g$b_BP)),
        stringsAsFactors = FALSE
      )
    }
  )
)

utils::write.csv(
  budget_all_methods,
  file.path(data_dir, "05_budget_all_methods.csv"),
  row.names = FALSE
)

write_tex_table(
  budget_all_methods %>%
    transmute(
      `k/n` = fmt_pct(k_fraction, 1),
      `MIS reach` = fmt_pct(MIS_reach, 1),
      `Normalized MIS boundary` = fmt_num(
        MIS_censored_median_normalized, 3
      ),
      `Cook reach` = fmt_pct(Cook_reach, 1),
      `DFBETAS reach` = fmt_pct(DFBETAS_reach, 1),
      `Leverage reach` = fmt_pct(Leverage_reach, 1),
      `RESET reach` = fmt_pct(RESET_reach, 1),
      `BP reach` = fmt_pct(BP_reach, 1)
    ),
  exhibit_path("tab05_budget"),
  caption = paste0(
    "Detection reach by deletion fraction. The pooled MIS boundary is ",
    "normalized by each scenario's maximum simulated severity before the ",
    "censoring-aware median is formed."
  ),
  label = "tab:05-budget-summary",
  resize_width = TABLE_WIDTH_WIDE,
  align = "lrrrrrrr",
  placement = "htbp"
)


budget_note <- paste0(
  "Boundaries are right-censored at the maximum simulated severity, so the ",
  "censored-sample median is reported; it is identified whenever fewer than ",
  "half of environments are censored. Reacher-only means and medians, together ",
  "with the censoring counts, are given in Table~\\ref{tab:05-budget-mean}."
)

# Body: scenario x k, censoring retained explicitly.
write_tex_table(
  budget_mis %>%
    mutate(
      Scenario = unname(scenario_labels[scenario]),
      scenario_order = factor(
        Scenario,
        levels = unname(scenario_labels[scenario_levels])
      )
    ) %>%
    arrange(scenario_order, k_fraction) %>%
    transmute(
      Scenario,
      `k/n` = fmt_pct(k_fraction, 1),
      `Reachers / 96` = paste0(reachers, " / ", total),
      Censored = censored,
      `Censoring-aware median` = fmt_censored_boundary(
        censored_median,
        median_identified,
        censor_limit,
        3
      )
    ),
  exhibit_path("tab_budget_by_scenario"),
  caption = "MIS deletion-budget sensitivity by scenario.",
  label = "tab:05-budget-by-scenario",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrr",
  placement = "htbp",
  note = budget_note,
  note_raw = TRUE,
  addlinespace_after = seq(4L, 28L, by = 4L)
)


# Appendix A5: same computation, with the conditional reacher mean/median.
write_tex_table(
  budget_mis %>%
    mutate(
      Scenario = unname(scenario_labels[scenario]),
      scenario_order = factor(
        Scenario,
        levels = unname(scenario_labels[scenario_levels])
      )
    ) %>%
    arrange(scenario_order, k_fraction) %>%
    transmute(
      Scenario,
      `k/n` = fmt_pct(k_fraction, 1),
      `Reachers / 96` = paste0(reachers, " / ", total),
      Censored = censored,
      `Reach rate` = fmt_pct(reach_rate, 1),
      `Censoring-aware median` = fmt_censored_boundary(
        censored_median,
        median_identified,
        censor_limit,
        3
      ),
      `Reacher mean` = fmt_num(reacher_mean, 3),
      `Reacher median` = fmt_num(reacher_median, 3)
    ),
  exhibit_path("tabA05_budget_mean"),
  caption = paste0(
    "Deletion-budget boundary audit. Reacher-only means and medians condition ",
    "on successful detection; censoring counts show the conditioning set ",
    "explicitly for every scenario and k."
  ),
  label = "tab:05-budget-mean",
  resize_width = TABLE_WIDTH_WIDE,
  align = "lrrrrrrr",
  placement = "htbp",
  addlinespace_after = seq(4L, 28L, by = 4L)
)


# Dynamic budget notes for drafting.
detectable_budget <- budget_mis %>%
  group_by(scenario) %>%
  filter(max(reach_rate, na.rm = TRUE) > 0) %>%
  mutate(
    best = max(reach_rate, na.rm = TRUE),
    primary_tied_best = any(
      near_num(k_fraction, PRIMARY_K) &
        abs(reach_rate - best) < 1e-12
    )
  ) %>%
  summarise(
    primary_tied_best = first(primary_tied_best),
    .groups = "drop"
  )

n_best_primary <- sum(detectable_budget$primary_tied_best)
n_detectable <- nrow(detectable_budget)

pooled_primary_row <- budget_all_methods %>%
  filter(near_num(k_fraction, PRIMARY_K))

budget_notes <- c(
  paste0(
    "Primary k/n = ", PRIMARY_K,
    " is best-or-tied on 80% reach in ",
    n_best_primary, " of ", n_detectable,
    " scenarios with any finite MIS boundary."
  ),
  paste0(
    "Pooled MIS reach at primary k/n: ",
    fmt_pct(pooled_primary_row$MIS_reach, 1),
    " (", pooled_primary_row$MIS_count, " / 768)."
  ),
  paste0(
    "Censoring-aware pooled normalized MIS median at primary k/n: ",
    fmt_num(pooled_primary_row$MIS_censored_median_normalized, 4),
    "."
  )
)

writeLines(
  budget_notes,
  file.path(diag_dir, "05_budget_notes.txt"),
  useBytes = TRUE
)


# ==============================================================================
# 8. SPECIFICITY & INVARIANCE
# ==============================================================================

# cell_id is exact and uniquely identifies the evaluation design cell.
# Do not join on severity_target or k_fraction: both are doubles and therefore
# inappropriate as identity keys.
wrong_power_for_spec <- W %>%
  filter(near_num(k_fraction, PRIMARY_K)) %>%
  select(
    cell_id,
    wrong_power = rej_MIS
  )

if (anyDuplicated(wrong_power_for_spec$cell_id)) {
  stop(
    "cell_id is not unique in the primary-k wrong-model detector summary."
  )
}

spec_primary <- SPEC %>%
  filter(near_num(k_fraction, PRIMARY_K)) %>%
  left_join(
    wrong_power_for_spec,
    by = "cell_id"
  )

if (any(is.na(spec_primary$wrong_power))) {
  warning(
    "Some 05c rows could not be matched to wrong-model MIS power. ",
    "Check cell identifiers if the specificity exhibit looks incomplete."
  )
}

spec_summary <- spec_primary %>%
  filter(severity_index > 1L) %>%
  group_by(scenario) %>%
  summarise(
    wrong_power = safe_mean(wrong_power),
    correct_fpr = safe_mean(power_correct_calibrated),
    mean_wrong_MIS = safe_mean(mean_MIS_wrong),
    mean_correct_MIS = safe_mean(mean_MIS_correct),
    .groups = "drop"
  ) %>%
  mutate(
    scenario_display = factor(
      unname(scenario_labels[scenario]),
      levels = rev(unname(scenario_labels[
        c(
          "nonlinear",
          "missing_interaction",
          "heterogeneous",
          "threshold",
          "ovb"
        )
      ]))
    )
  )

utils::write.csv(
  spec_summary,
  file.path(data_dir, "05_specificity_summary_primary_k.csv"),
  row.names = FALSE
)

spec_plot_data <- spec_summary %>%
  select(scenario_display, wrong_power, correct_fpr) %>%
  pivot_longer(
    cols = c(wrong_power, correct_fpr),
    names_to = "state",
    values_to = "rejection_rate"
  ) %>%
  mutate(
    state = factor(
      state,
      levels = c("wrong_power", "correct_fpr"),
      labels = c("Wrong model", "Correct model")
    )
  )

fig_specificity <- ggplot(
  spec_plot_data,
  aes(
    x = rejection_rate,
    y = scenario_display,
    fill = state
  )
) +
  geom_vline(
    xintercept = 0.05,
    colour = COL_GREY_DARK,
    linetype = "dashed",
    linewidth = 0.60
  ) +
  geom_col(
    position = position_dodge2(width = 0.78, preserve = "single"),
    width = 0.62,
    alpha = 0.92
  ) +
  scale_fill_manual(
    values = c(
      "Wrong model" = COL_GREY,
      "Correct model" = COL_MIS
    ),
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Rejection rate",
    y = NULL
  ) +
  theme_85() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )

save_plot(
  fig_specificity,
  exhibit_path("fig08_specificity"),
  width = 7.8,
  height = 4.9
)


write_tex_table(
  spec_summary %>%
    transmute(
      Scenario = as.character(scenario_display),
      `Wrong-model rejection` = fmt_pct(wrong_power, 1),
      `Correct-model rejection` = fmt_pct(correct_fpr, 2)
    ),
  exhibit_path("tab06_specificity"),
  caption = paste0(
    "Wrong-model rejection probability and correctly calibrated false-positive ",
    "rate at k/n = 2.5 percent for the five scenarios with an explicit ",
    "corrected-model specification in 05c."
  ),
  label = "tab:05-specificity",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrr",
  placement = "htbp"
)


# ------------------------------------------------------------------------------
# Invariance / blind-set diagnostics
# ------------------------------------------------------------------------------

endo <- P %>%
  filter(scenario == "endogeneity") %>%
  arrange(seed_group_id, severity_index)

endo_dev <- endo %>%
  group_by(seed_group_id) %>%
  summarise(
    max_abs_deviation = max(
      abs(mean_MIS - first(mean_MIS)),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

endo_positive <- endo %>%
  filter(severity_index > 1L)

endo_max <- endo %>%
  group_by(env_id) %>%
  filter(severity_index == max(severity_index)) %>%
  ungroup()

ovb <- P %>%
  filter(scenario == "ovb")

ovb_positive <- ovb %>%
  filter(severity_index > 1L)

ovb_max <- ovb %>%
  group_by(env_id) %>%
  filter(severity_index == max(severity_index)) %>%
  ungroup()

correct_spread <- spec_primary %>%
  group_by(env_id, scenario) %>%
  summarise(
    statistic_spread = max(mean_MIS_correct, na.rm = TRUE) -
      min(mean_MIS_correct, na.rm = TRUE),
    correct_fpr = safe_mean(power_correct_calibrated),
    .groups = "drop"
  )

blind_table_numeric <- data.frame(
  case = c(
    "Linear endogeneity",
    "Omitted variable",
    "Correctly respecified models"
  ),
  interpretation = c(
    "Exact fitted-span / residual-scale invariance",
    "Near-blind case: weak residual-geometry perturbation",
    "Specificity control after adding the missing structure"
  ),
  mean_positive_rejection = c(
    safe_mean(endo_positive$rej_MIS),
    safe_mean(ovb_positive$rej_MIS),
    safe_mean(correct_spread$correct_fpr)
  ),
  max_severity_rejection = c(
    safe_mean(endo_max$rej_MIS),
    safe_mean(ovb_max$rej_MIS),
    NA_real_
  ),
  max_statistic_deviation = c(
    max(endo_dev$max_abs_deviation, na.rm = TRUE),
    NA_real_,
    max(correct_spread$statistic_spread, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  blind_table_numeric,
  file.path(data_dir, "05_blind_invariance_summary.csv"),
  row.names = FALSE
)

write_tex_table(
  blind_table_numeric %>%
    transmute(
      Case = case,
      Interpretation = interpretation,
      `Mean rejection` = fmt_pct(mean_positive_rejection, 2),
      `Max-severity rejection` = fmt_pct(max_severity_rejection, 2),
      `Max statistic deviation` = fmt_sci(max_statistic_deviation, 2)
    ),
  exhibit_path("tab07_invariance"),
  caption = "Blind-set and invariance diagnostics at k/n = 2.5 percent.",
  label = "tab:05-invariance-blind-set",
  resize_width = TABLE_WIDTH_WIDE,
  align = "llrrr",
  placement = "htbp",
  note = paste0(
    "A machine-scale maximum statistic deviation under linear endogeneity ",
    "supports the fitted-span invariance result rather than merely indicating ",
    "low power."
  )
)


# Algebraic argument as a ready-to-input LaTeX text fragment.
invariance_tex <- c(
  "% Auto-generated by scripts/85_output_misspecification.R",
  "\\paragraph{Why the fitted-span blind case is invariant.}",
  paste0(
    "Suppose the severity-indexed data can be written as ",
    "\\(y_{\\lambda}=X\\beta+a_{\\lambda}X\\gamma+b_{\\lambda}u\\). ",
    "For the full sample, \\(\\widehat\\beta_{\\lambda}=\\beta+",
    "a_{\\lambda}\\gamma+b_{\\lambda}(X'X)^{-1}X'u\\). ",
    "After deleting any fixed set \\(S\\), the fitted coefficient has the same ",
    "additive \\(a_{\\lambda}\\gamma\\) shift, so the full-versus-deleted ",
    "coefficient difference is multiplied only by \\(b_{\\lambda}\\). "
  ),
  paste0(
    "The OLS standard error is also multiplied by \\(|b_{\\lambda}|\\). ",
    "Therefore the standardized deletion statistic satisfies"
  ),
  "\\[",
  "\\frac{\\left|\\widehat\\beta_{\\lambda}-\\widehat\\beta_{\\lambda,-S}\\right|}",
  "{\\operatorname{SE}(\\widehat\\beta_{\\lambda})}",
  "=",
  "\\frac{\\left|\\widehat\\beta_{0}-\\widehat\\beta_{0,-S}\\right|}",
  "{\\operatorname{SE}(\\widehat\\beta_{0})}.",
  "\\]",
  paste0(
    "Thus a departure that changes only the fitted-span coefficient and an ",
    "overall residual scale is invisible to the standardized MIS statistic. ",
    "Because this equality holds for every admissible deletion set \\(S\\), ",
    "maximizing the standardized statistic over \\(S\\) preserves the same ",
    "invariance. The common-random-number simulation checks this proposition ",
    "numerically."
  )
)

writeLines(
  invariance_tex,
  exhibit_path("text07_invariance_argument"),
  useBytes = TRUE
)


# ==============================================================================
# 9. LOCALISATION
# ==============================================================================

localisation <- P %>%
  filter(
    severity_index > 1L,
    scenario %in% c("threshold", "heterogeneous"),
    is.finite(mean_precision),
    is.finite(mean_recall),
    is.finite(mean_lift),
    is.finite(mean_affected)
  ) %>%
  mutate(
    # 05a uses round(k_fraction * n), so the realised selected fraction is
    # slightly below 2.5% for n = 500 and 2500 under R's rounding rule.
    selected_fraction = pmax(1, round(k_fraction * n)) / n,
    recall_ceiling_row = pmin(
      1,
      selected_fraction / mean_affected
    ),
    recall_efficiency_row = mean_recall / recall_ceiling_row,
    recall_precision_ratio_row = mean_recall / mean_precision
  ) %>%
  group_by(scenario) %>%
  summarise(
    affected_fraction = safe_mean(mean_affected),
    selected_fraction = safe_mean(selected_fraction),
    precision = safe_mean(mean_precision),
    recall = safe_mean(mean_recall),
    lift = safe_mean(mean_lift),
    recall_ceiling = safe_mean(recall_ceiling_row),
    recall_efficiency = safe_mean(recall_efficiency_row),
    recall_precision_ratio = safe_mean(recall_precision_ratio_row),
    .groups = "drop"
  ) %>%
  mutate(
    scenario_display = factor(
      unname(scenario_labels[scenario]),
      levels = c("Heterogeneous slope", "Threshold")
    )
  )

utils::write.csv(
  localisation,
  file.path(data_dir, "05_localisation_primary_k.csv"),
  row.names = FALSE
)

write_tex_table(
  localisation %>%
    transmute(
      Scenario = as.character(scenario_display),
      Precision = fmt_num(precision, 3),
      Lift = fmt_num(lift, 2)
    ),
  exhibit_path("tab08_localisation"),
  caption = "MIS localisation at k/n = 2.5 percent.",
  label = "tab:05-localisation",
  resize_width = TABLE_WIDTH_COMPACT,
  align = "lrr",
  placement = "htbp",
  note = paste0(
    "The affected fraction is 25 percent. Lift is precision relative to this ",
    "chance benchmark. Because recall efficiency reduces algebraically to ",
    "precision under the fixed affected-set design, it is not reported separately."
  )
)


loc_plot_data <- bind_rows(
  localisation %>%
    transmute(
      scenario_display,
      metric = "Precision",
      value = precision,
      reference = affected_fraction
    ),
  localisation %>%
    transmute(
      scenario_display,
      metric = "Lift",
      value = lift,
      reference = 1
    )
)

fig_localisation <- ggplot(
  loc_plot_data,
  aes(
    x = value,
    y = scenario_display
  )
) +
  geom_segment(
    aes(
      x = reference,
      xend = value,
      yend = scenario_display
    ),
    linewidth = 1.4,
    colour = COL_BLUE_PALE
  ) +
  geom_point(
    size = 3.0,
    shape = 21,
    fill = COL_MIS,
    colour = COL_BLUE_DARK,
    stroke = 0.55
  ) +
  geom_point(
    aes(x = reference),
    size = 2.4,
    shape = 23,
    fill = COL_NEUTRAL,
    colour = COL_GREY_DARK,
    stroke = 0.5
  ) +
  facet_wrap(
    ~ metric,
    nrow = 1,
    scales = "free_x"
  ) +
  labs(
    x = "Observed value; diamond marks chance/reference level",
    y = NULL
  ) +
  theme_85() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# fig_localisation is intentionally not saved alone. It is combined with the
# deletion-arm ladder in Section 10; if frozen 05 is unavailable, Section 10
# falls back to saving this panel by itself.


# ==============================================================================
# 10. DELETION-ARM LADDER — FROZEN 05, TARGETING NOT CLEANING
# ==============================================================================

if (HAS_FROZEN_05) {
  
  F05 <- readRDS(input_05_frozen_summary)
  
  deletion_ladder <- F05 %>%
    filter(
      model_state == "wrong",
      estimator == "OLS",
      # Frozen 05 is 1-indexed: severity_level == 1 is the null cell.
      severity_level > 1L,
      diagnostic %in% c("none", "MIS", "Cook", "DFBETAS", "Leverage")
    ) %>%
    mutate(
      arm = case_when(
        diagnostic == "none" ~ "Full OLS",
        diagnostic == "MIS" ~ "MIS deletion",
        diagnostic == "Cook" ~ "Cook deletion",
        diagnostic == "DFBETAS" ~ "DFBETAS deletion",
        diagnostic == "Leverage" ~ "Leverage deletion",
        TRUE ~ diagnostic
      )
    ) %>%
    group_by(x_type, k_fraction, arm) %>%
    summarise(
      mean_cell_RMSE = safe_mean(rmse),
      cells = n(),
      .groups = "drop"
    ) %>%
    mutate(
      x_display = factor(
        unname(x_labels[x_type]),
        levels = unname(x_labels[
          c("normal", "mixed_normal", "contaminated")
        ])
      ),
      arm = factor(
        arm,
        levels = c(
          "Full OLS",
          "MIS deletion",
          "Cook deletion",
          "DFBETAS deletion",
          "Leverage deletion"
        )
      )
    )
  
  utils::write.csv(
    deletion_ladder,
    file.path(data_dir, "05_frozen_deletion_arm_ladder_by_x.csv"),
    row.names = FALSE
  )
  
  frozen_rmse_mean_median <- F05 %>%
    filter(
      model_state == "wrong",
      estimator == "OLS",
      severity_level > 1L,
      diagnostic %in% c("none", "MIS", "Cook", "DFBETAS", "Leverage")
    ) %>%
    mutate(
      arm = case_when(
        diagnostic == "none" ~ "Full OLS",
        diagnostic == "MIS" ~ "MIS deletion",
        diagnostic == "Cook" ~ "Cook deletion",
        diagnostic == "DFBETAS" ~ "DFBETAS deletion",
        diagnostic == "Leverage" ~ "Leverage deletion",
        TRUE ~ diagnostic
      )
    ) %>%
    group_by(x_type, k_fraction, arm) %>%
    summarise(
      mean_cell_RMSE = safe_mean(rmse),
      median_cell_RMSE = safe_median(rmse),
      cells = n(),
      .groups = "drop"
    )
  
  utils::write.csv(
    frozen_rmse_mean_median,
    file.path(data_dir, "05_frozen_rmse_mean_median.csv"),
    row.names = FALSE
  )
  
  write_tex_table(
    frozen_rmse_mean_median %>%
      mutate(
        x_order = factor(
          x_type,
          levels = c(
            "normal",
            "mixed_normal",
            "contaminated"
          )
        ),
        arm_order = factor(
          arm,
          levels = c(
            "Full OLS",
            "MIS deletion",
            "Cook deletion",
            "DFBETAS deletion",
            "Leverage deletion"
          )
        )
      ) %>%
      arrange(
        x_order,
        k_fraction,
        arm_order
      ) %>%
      transmute(
        `X distribution` = unname(x_labels[x_type]),
        `k/n` = fmt_pct(k_fraction, 1),
        Arm = arm,
        `Mean cell RMSE` = fmt_num(mean_cell_RMSE, 3),
        `Median cell RMSE` = fmt_num(median_cell_RMSE, 3)
      ),
    exhibit_path("tabA06_frozen_rmse_mean"),
    caption = paste0(
      "Frozen-05 deletion-arm RMSE sensitivity. Arithmetic means and medians ",
      "of the cell-level RMSE are shown side by side; the severity-0 null cell ",
      "is excluded and X distributions remain separated."
    ),
    label = "tab:05-frozen-rmse-mean-median",
    resize_width = TABLE_WIDTH_WIDE,
    align = "lllrr",
    placement = "htbp",
    addlinespace_after = seq(5L, 55L, by = 5L)
  )
  
  arm_colours <- c(
    "Full OLS" = COL_GREY_DARK,
    "MIS deletion" = COL_MIS,
    "Cook deletion" = COL_COOK,
    "DFBETAS deletion" = COL_DFBETAS,
    "Leverage deletion" = COL_LEVERAGE
  )
  
  fig_ladder <- ggplot(
    deletion_ladder,
    aes(
      x = k_fraction,
      y = mean_cell_RMSE,
      colour = arm,
      linetype = arm,
      group = arm
    )
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.0) +
    facet_wrap(
      ~ x_display,
      nrow = 1,
      scales = "free_y"
    ) +
    scale_x_continuous(
      breaks = sort(unique(deletion_ladder$k_fraction)),
      labels = function(x) scales::percent(x, accuracy = 0.1)
    ) +
    scale_y_continuous(
      transform = scales::pseudo_log_trans(base = 10, sigma = 0.20)
    ) +
    scale_colour_manual(values = arm_colours) +
    scale_linetype_manual(
      values = c(
        "Full OLS" = "dashed",
        "MIS deletion" = "solid",
        "Cook deletion" = "22",
        "DFBETAS deletion" = "42",
        "Leverage deletion" = "13"
      )
    ) +
    labs(
      x = "Deletion fraction k/n",
      y = "Mean cell-level RMSE (pseudo-log scale)",
      colour = NULL,
      linetype = NULL
    ) +
    theme_85(base_size = 9.6) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 8.5)
    )
  
  save_stacked_plots(
    top_plot = fig_localisation,
    bottom_plot = fig_ladder,
    pdf_path = exhibit_path("fig09_targeting_localisation"),
    width = 9.4,
    height = 8.4,
    top_share = 0.36
  )
  
} else {
  F05 <- NULL
  
  # The formal localisation panel remains useful even if the frozen-05
  # estimation summary is unavailable.
  save_plot(
    fig_localisation,
    exhibit_path("fig09_targeting_localisation"),
    width = 8.0,
    height = 3.5
  )
}


# ==============================================================================
# 11. ENVIRONMENT HETEROGENEITY
# ==============================================================================

power_summary_by <- function(group_var) {
  PW %>%
    group_by(
      across(all_of(group_var))
    ) %>%
    summarise(
      MIS = safe_mean(rej_MIS),
      Cook = safe_mean(rej_Cook),
      DFBETAS = safe_mean(rej_DFBETAS),
      Leverage = safe_mean(rej_Leverage),
      RESET = safe_mean(rej_RESET_cal),
      BP = safe_mean(rej_BP_cal),
      .groups = "drop"
    )
}

env_error <- power_summary_by("error_type")
env_x <- power_summary_by("x_type")
env_n <- power_summary_by("n")

utils::write.csv(
  env_error,
  file.path(data_dir, "05_environment_power_by_error.csv"),
  row.names = FALSE
)

utils::write.csv(
  env_x,
  file.path(data_dir, "05_environment_power_by_x.csv"),
  row.names = FALSE
)

utils::write.csv(
  env_n,
  file.path(data_dir, "05_environment_power_by_n.csv"),
  row.names = FALSE
)


write_tex_table(
  env_error %>%
    transmute(
      `Error distribution` = unname(error_labels[error_type]),
      MIS = fmt_pct(MIS, 1),
      `Cook's D` = fmt_pct(Cook, 1),
      DFBETAS = fmt_pct(DFBETAS, 1),
      Leverage = fmt_pct(Leverage, 1),
      RESET = fmt_pct(RESET, 1),
      BP = fmt_pct(BP, 1)
    ),
  exhibit_path("tabS05_power_error"),
  caption = paste0(
    "Mean detection power across positive-severity simulation cells by error ",
    "distribution at k/n = 2.5 percent; scenario-specific severity grids are ",
    "retained."
  ),label = "tab:05-power-by-error",
  resize_width = TABLE_WIDTH_WIDE,
  align = "lrrrrrr",
  placement = "htbp"
)


write_tex_table(
  env_x %>%
    mutate(
      x_order = factor(
        x_type,
        levels = c("normal", "mixed_normal", "contaminated")
      )
    ) %>%
    arrange(x_order) %>%
    transmute(
      `X distribution` = unname(x_labels[x_type]),
      MIS = fmt_pct(MIS, 1),
      `Cook's D` = fmt_pct(Cook, 1),
      DFBETAS = fmt_pct(DFBETAS, 1),
      Leverage = fmt_pct(Leverage, 1),
      RESET = fmt_pct(RESET, 1),
      BP = fmt_pct(BP, 1)
    ),
  exhibit_path("tabS06_power_x"),
  caption = paste0(
    "Mean detection power across positive-severity simulation cells by X ",
    "distribution at k/n = 2.5 percent; scenario-specific severity grids are ",
    "retained."
  ),
  label = "tab:05-power-by-x",
  resize_width = TABLE_WIDTH_WIDE,
  align = "lrrrrrr",
  placement = "htbp"
)


write_tex_table(
  env_n %>%
    transmute(
      n = n,
      MIS = fmt_pct(MIS, 1),
      `Cook's D` = fmt_pct(Cook, 1),
      DFBETAS = fmt_pct(DFBETAS, 1),
      Leverage = fmt_pct(Leverage, 1),
      RESET = fmt_pct(RESET, 1),
      BP = fmt_pct(BP, 1)
    ),
  exhibit_path("tabS07_power_n"),
  caption = paste0(
    "Mean detection power across positive-severity simulation cells by sample ",
    "size at k/n = 2.5 percent; scenario-specific severity grids are retained."
  ),
  label = "tab:05-power-by-n",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrrrr",
  placement = "htbp"
)


# Across-environment power distribution at each scenario's maximum severity.
max_severity_rows <- P %>%
  group_by(scenario, env_id) %>%
  filter(severity_index == max(severity_index)) %>%
  ungroup()

worst_case <- max_severity_rows %>%
  group_by(scenario) %>%
  summarise(
    min_power = min(rej_MIS, na.rm = TRUE),
    p10 = safe_quantile(rej_MIS, 0.10),
    median = safe_median(rej_MIS),
    max_power = max(rej_MIS, na.rm = TRUE),
    .groups = "drop"
  )

utils::write.csv(
  worst_case,
  file.path(data_dir, "05_worst_case_max_severity.csv"),
  row.names = FALSE
)

write_tex_table(
  worst_case %>%
    transmute(
      Scenario = unname(scenario_labels[scenario]),
      Minimum = fmt_pct(min_power, 1),
      `10th percentile` = fmt_pct(p10, 1),
      Median = fmt_pct(median, 1),
      Maximum = fmt_pct(max_power, 1)
    ),
  exhibit_path("tabS08_worst_environment"),
  caption = paste0(
    "Distribution of MIS power across the 96 environments at each scenario's ",
    "maximum simulated severity, k/n = 2.5 percent."
  ),
  label = "tab:05-worst-case-environment",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrr",
  placement = "htbp"
)


# Soft heatmap for environment heterogeneity.
env_heat <- PW %>%
  group_by(n, x_type, error_type) %>%
  summarise(
    MIS_power = safe_mean(rej_MIS),
    .groups = "drop"
  ) %>%
  mutate(
    x_display = factor(
      unname(x_labels[x_type]),
      levels = unname(x_labels[
        c("normal", "mixed_normal", "contaminated")
      ])
    ),
    error_display = factor(
      unname(error_labels[error_type]),
      levels = rev(unname(error_labels[
        c(
          "normal", "mixed_normal", "skewed_t", "golm",
          "beta_logistic", "gpd", "contaminated", "pareto"
        )
      ]))
    )
  )

utils::write.csv(
  env_heat,
  file.path(data_dir, "05_environment_heatmap_data.csv"),
  row.names = FALSE
)

fig_env <- ggplot(
  env_heat,
  aes(
    x = x_display,
    y = error_display,
    fill = MIS_power
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.75
  ) +
  geom_text(
    aes(
      label = fmt_pct(MIS_power, 0)
    ),
    size = 2.6,
    colour = COL_GREY_DARK
  ) +
  facet_wrap(
    ~ n,
    nrow = 1
  ) +
  scale_fill_gradient(
    low = COL_NEUTRAL,
    high = COL_MIS,
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = scales::percent_format(accuracy = 1),
    guide = guide_colourbar(
      barwidth = grid::unit(7, "cm"),
      barheight = grid::unit(0.45, "cm")
    )
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "MIS power"
  ) +
  theme_85_heatmap(base_size = 9.2) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      vjust = 1
    ),
    strip.text = element_text(size = 9.5, face = "bold")
  )

save_plot(
  fig_env,
  exhibit_path("figS02_environment"),
  width = 11.0,
  height = 5.0
)


# ==============================================================================
# 12. FROZEN 05 APPENDIX — ESTIMATION ONLY
# ==============================================================================

if (HAS_FROZEN_05) {
  
  # ---------------------------------------------------------------------------
  # 12.1 Correct-vs-wrong MM: bias, RMSE, coverage, SE ratio, predicted coverage
  #
  # Full MM is repeated for each k in frozen 05. Keep one k only to avoid
  # duplicate weighting. This DOES NOT give k any substantive role for full MM.
  # ---------------------------------------------------------------------------
  
  # Frozen 05 is 1-indexed: severity_level == 1 is severity 0.
  #
  # Wrong-state estimation summaries exclude that null cell (> 1L).
  # Correct-state rows deliberately retain the severity-0 control (> 0L):
  # the correctly specified model is used here as a calibration/control branch,
  # and its near-nominal coverage across the full branch is itself informative.
  mm05_wrong <- F05 %>%
    filter(
      near_num(k_fraction, PRIMARY_K),
      diagnostic == "none",
      estimator == "MM",
      model_state == "wrong",
      severity_level > 1L
    )
  
  mm05_correct <- F05 %>%
    filter(
      near_num(k_fraction, PRIMARY_K),
      diagnostic == "none",
      estimator == "MM",
      model_state == "correct",
      severity_level > 0L
    )
  
  mm05 <- bind_rows(mm05_wrong, mm05_correct) %>%
    mutate(
      sampling_sd_proxy = sqrt(pmax(rmse^2 - mean_bias^2, 0)),
      se_ratio = ifelse(
        is.finite(sampling_sd_proxy) &
          sampling_sd_proxy > sqrt(.Machine$double.eps),
        mean_se / sampling_sd_proxy,
        NA_real_
      ),
      predicted_coverage = ifelse(
        is.finite(sampling_sd_proxy) &
          sampling_sd_proxy > sqrt(.Machine$double.eps) &
          is.finite(mean_se) &
          is.finite(mean_bias),
        stats::pnorm(
          (1.96 * mean_se - mean_bias) / sampling_sd_proxy
        ) -
          stats::pnorm(
            (-1.96 * mean_se - mean_bias) / sampling_sd_proxy
          ),
        NA_real_
      )
    )
  
  # Only scenarios where both states are actually available.
  valid_state_scenarios <- mm05 %>%
    distinct(scenario, model_state) %>%
    count(scenario, name = "n_states") %>%
    filter(n_states >= 2L) %>%
    pull(scenario)
  
  mm05_compare <- mm05 %>%
    filter(scenario %in% valid_state_scenarios) %>%
    group_by(x_type, scenario, model_state) %>%
    summarise(
      mean_abs_bias = safe_mean(mean_abs_bias),
      rmse = safe_mean(rmse),
      coverage = safe_mean(coverage),
      se_ratio = safe_mean(se_ratio),
      predicted_coverage = safe_mean(predicted_coverage),
      cells = n(),
      .groups = "drop"
    )
  
  utils::write.csv(
    mm05_compare,
    file.path(data_dir, "05_frozen_mm_correct_vs_wrong.csv"),
    row.names = FALSE
  )
  
  write_tex_table(
    mm05_compare %>%
      mutate(
        x_order = factor(
          x_type,
          levels = c("normal", "mixed_normal", "contaminated")
        )
      )%>%
      transmute(
        `X distribution` = unname(x_labels[x_type]),
        Scenario = ifelse(
          scenario %in% names(frozen_scenario_labels),
          unname(frozen_scenario_labels[scenario]),
          scenario
        ),
        State = ifelse(
          model_state == "correct",
          "Correct",
          "Wrong"
        ),
        `Abs. bias` = fmt_num(mean_abs_bias, 3),
        RMSE = fmt_num(rmse, 3),
        Coverage = fmt_pct(coverage, 1),
        `SE ratio` = fmt_num(se_ratio, 3),
        `Predicted coverage` = fmt_pct(predicted_coverage, 1)
      ) %>%
      arrange(`X distribution`, Scenario, State),
    exhibit_path("tabA01_mm_correct_wrong"),
    caption = paste0(
      "Frozen-05 full-MM estimation under wrong and correctly respecified ",
      "models, reported separately by X distribution."
    ),
    label = "tab:05-frozen-mm-correct-wrong",
    resize_width = TABLE_WIDTH_WIDE,
    align = "lllrrrrr",
    placement = "htbp",
    note = paste0(
      "SE ratio is mean reported SE divided by the sampling-SD proxy ",
      "sqrt(RMSE squared minus mean bias squared); predicted coverage is the ",
      "corresponding normal approximation."
    ),
    addlinespace_after = c(12L, 24L)
  )

  # ---------------------------------------------------------------------------
  # 12.2 Full OLS vs full MM, by x_type and scenario
  # ---------------------------------------------------------------------------
  
  base05 <- F05 %>%
    filter(
      near_num(k_fraction, PRIMARY_K),
      model_state == "wrong",
      diagnostic == "none",
      estimator %in% c("OLS", "MM"),
      # Exclude frozen-05 null cell: severity_level == 1 is severity 0.
      severity_level > 1L
    ) %>%
    group_by(x_type, scenario, estimator) %>%
    summarise(
      mean_abs_bias = safe_mean(mean_abs_bias),
      rmse = safe_mean(rmse),
      coverage = safe_mean(coverage),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = estimator,
      values_from = c(mean_abs_bias, rmse, coverage),
      names_glue = "{.value}_{estimator}"
    ) %>%
    mutate(
      MM_lower_RMSE = rmse_MM < rmse_OLS,
      MM_closer_coverage = abs(coverage_MM - 0.95) <
        abs(coverage_OLS - 0.95)
    )
  
  utils::write.csv(
    base05,
    file.path(data_dir, "05_frozen_ols_vs_mm_by_x_scenario.csv"),
    row.names = FALSE
  )
  
  write_tex_table(
    base05 %>%
      mutate(
        x_order = factor(
          x_type,
          levels = c("normal", "mixed_normal", "contaminated")
        )
      ) %>%
      arrange(x_order, scenario) %>%
      transmute(
        `X distribution` = unname(x_labels[x_type]),
        Scenario = ifelse(
          scenario %in% names(frozen_scenario_labels),
          unname(frozen_scenario_labels[scenario]),
          scenario
        ),
        `OLS RMSE` = fmt_num(rmse_OLS, 3),
        `MM RMSE` = fmt_num(rmse_MM, 3),
        `OLS coverage` = fmt_pct(coverage_OLS, 1),
        `MM coverage` = fmt_pct(coverage_MM, 1),
        `MM closer to 95 percent` = ifelse(
          MM_closer_coverage, "Yes", "No"
        )
      ),
    exhibit_path("tabA02_ols_mm"),
    caption = paste0(
      "Frozen-05 estimation comparison of full OLS and full MM. ",
      "Results are reported within X-distribution strata; no comparison of ",
      "frozen-05 severity values with the formal 05a scale is made."
    ),
    label = "tab:05-frozen-ols-vs-mm",
    resize_width = TABLE_WIDTH_WIDE,
    align = "llrrrrl",
    placement = "htbp",
    addlinespace_after = c(8L, 16L)
  )
  
  
  mm_win_counts <- base05 %>%
    group_by(x_type) %>%
    summarise(
      scenarios = n(),
      MM_RMSE_wins = sum(MM_lower_RMSE, na.rm = TRUE),
      MM_coverage_wins = sum(MM_closer_coverage, na.rm = TRUE),
      .groups = "drop"
    )
  
  utils::write.csv(
    mm_win_counts,
    file.path(data_dir, "05_frozen_mm_win_counts_by_x.csv"),
    row.names = FALSE
  )
  
  write_tex_table(
    mm_win_counts %>%
      transmute(
        `X distribution` = unname(x_labels[x_type]),
        `MM RMSE wins` = paste0(MM_RMSE_wins, "/", scenarios),
        `MM coverage wins` = paste0(MM_coverage_wins, "/", scenarios)
      ),
    exhibit_path("tabA03_mm_wins"),
    caption = paste0(
      "Frozen-05 count of scenarios in which full MM has lower RMSE or coverage ",
      "closer to 95 percent than full OLS, reported separately by X distribution."
    ),
    label = "tab:05-frozen-mm-win-counts",
    resize_width = TABLE_WIDTH_COMPACT,
    align = "lrr",
    placement = "htbp"
  )
}


# ==============================================================================
# 13. CORRECT-MODEL CUTOFF AUDIT
# ==============================================================================

cutoff_compare <- CORRECT_CUT %>%
  filter(near_num(k_fraction, PRIMARY_K)) %>%
  group_by(correct_class) %>%
  summarise(
    median_correct_cutoff = safe_median(cut_MIS),
    p10_correct_cutoff = safe_quantile(cut_MIS, 0.10),
    p90_correct_cutoff = safe_quantile(cut_MIS, 0.90),
    finite_rate = safe_mean(finite_MIS_rate),
    .groups = "drop"
  )

utils::write.csv(
  cutoff_compare,
  file.path(data_dir, "05_correct_model_cutoff_audit.csv"),
  row.names = FALSE
)

write_tex_table(
  cutoff_compare %>%
    transmute(
      `Correct model class` = ifelse(
        correct_class %in% names(scenario_labels),
        unname(scenario_labels[correct_class]),
        correct_class
      ),
      `Median cutoff` = fmt_num(median_correct_cutoff, 3),
      `10th percentile` = fmt_num(p10_correct_cutoff, 3),
      `90th percentile` = fmt_num(p90_correct_cutoff, 3),
      `Finite rate` = fmt_pct(finite_rate, 1)
    ),
  exhibit_path("tabS09_correct_cutoff"),
  caption = "Correct-model MIS null cutoffs from the independent 05c calibration at k/n = 2.5 percent.",
  label = "tab:05-correct-cutoff-audit",
  resize_width = TABLE_WIDTH_MEDIUM,
  align = "lrrrr",
  placement = "htbp"
)


# ==============================================================================
# 15. FINAL CONSOLE SUMMARY
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(diag_dir, "05_output_session_info.txt")
)

utils::write.csv(
  exhibit_registry,
  file.path(diag_dir, "05_exhibit_registry.csv"),
  row.names = FALSE
)


# Manifest is written last.
all_outputs <- list.files(
  output_root,
  recursive = TRUE,
  full.names = TRUE,
  all.files = FALSE
)

relative_to_output_root <- function(path) {
  root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)
  p <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")
  
  if (startsWith(p, prefix)) {
    substring(p, nchar(prefix) + 1L)
  } else {
    basename(p)
  }
}

manifest <- data.frame(
  relative_path = vapply(
    all_outputs,
    relative_to_output_root,
    character(1)
  ),
  bytes = file.info(all_outputs)$size,
  modified = as.character(file.info(all_outputs)$mtime),
  stringsAsFactors = FALSE
) %>%
  arrange(relative_path)

utils::write.csv(
  manifest,
  file.path(diag_dir, "05_output_manifest.csv"),
  row.names = FALSE
)


# ==============================================================================
# 14. FINAL CONSOLE SUMMARY
# ==============================================================================

cat("\n", strrep("=", 88), "\n", sep = "")
cat("SCRIPT 85 — MISSPECIFICATION OUTPUT COMPLETE\n")
cat(strrep("=", 88), "\n", sep = "")

cat("Output root: ", output_root, "\n", sep = "")
cat("Primary k/n: ", PRIMARY_K, "\n", sep = "")

cat("\nPrimary 80% reach:\n")
print(
  reach_counts %>%
    arrange(desc(reached)) %>%
    mutate(
      method = as.character(method),
      reach = paste0(reached, " / ", total, " (", fmt_pct(reach_rate, 1), ")")
    ) %>%
    select(method, reach),
  row.names = FALSE
)

cat("\nBoundary log-log slopes:\n")
print(
  boundary_loglog_slopes %>%
    transmute(
      scenario = scenario_display,
      n_points = n_points,
      slope = round(slope, 3),
      r_squared = round(r_squared, 3)
    ),
  row.names = FALSE
)

cat("\nBudget summary:\n")
print(
  budget_all_methods %>%
    transmute(
      k = fmt_pct(k_fraction, 1),
      MIS = fmt_pct(MIS_reach, 1),
      Cook = fmt_pct(Cook_reach, 1),
      DFBETAS = fmt_pct(DFBETAS_reach, 1),
      Leverage = fmt_pct(Leverage_reach, 1),
      RESET = fmt_pct(RESET_reach, 1),
      BP = fmt_pct(BP_reach, 1)
    ),
  row.names = FALSE
)

cat("\nSpecificity summary:\n")
print(
  spec_summary %>%
    transmute(
      scenario = as.character(scenario_display),
      wrong = fmt_pct(wrong_power, 1),
      correct = fmt_pct(correct_fpr, 2)
    ),
  row.names = FALSE
)

cat("\nLocalisation summary:\n")
print(
  localisation %>%
    transmute(
      scenario = as.character(scenario_display),
      precision = fmt_num(precision, 3),
      lift = fmt_num(lift, 2),
      recall_ceiling = fmt_num(recall_ceiling, 3),
      recall_efficiency = fmt_pct(recall_efficiency, 1)
    ),
  row.names = FALSE
)

cat("\nManifest:\n  ",
    file.path(diag_dir, "05_output_manifest.csv"), "\n", sep = "")