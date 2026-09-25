# ==============================================================================
# File: scripts/90_output_application.R
#
# Purpose:
#   Generate the cross-study empirical-application outputs for Section 5.
#
# Inputs:
#   application/<study>/audit_path.csv
#   application/<study>/audit_summary.csv
#   application/<study>/audit_baseline.csv
#
# Main outputs:
#   application/00_cross_study/
#     figures/main/
#       05_fig1_sensitivity_thresholds.pdf
#       05_fig2_influential_set_persistence.pdf
#
#     tables/main/
#       05_tab1_cross_study_sensitivity.tex
#       05_tab1_cross_study_sensitivity.csv
#
#     data/
#       05_cross_study_paths.csv
#       05_cross_study_thresholds.csv
#       05_cross_zero_mis_ids.csv
#       05_cross_zero_mis_ids_long.csv
#
#     diagnostics/
#       05_ratio_range_diagnostics.csv
#       05_figure3_status.txt
#       session_info.txt
#       manifest.txt
#
# Definitions:
#
#   R_k = beta_after / beta_original
#
#   50% attenuation:
#       R_k <= 0.5
#
#   50% amplification:
#       R_k >= 1.5
#
#   zero crossing / sign reversal:
#       R_k <= 0
#
# Interpretation:
#   - R_k = 1   : unchanged from baseline
#   - R_k = 0.5 : coefficient magnitude reduced by at least 50%
#   - R_k = 1.5 : coefficient magnitude increased by at least 50%
#                  in the original direction
#   - R_k = 0   : zero crossing
#   - R_k < 0   : sign reversal
#
# Figure 1:
#   Three threshold markers per study:
#     1. 50% attenuation
#     2. 50% amplification
#     3. zero crossing / sign reversal
#
#   Hollow markers at the right boundary indicate that the corresponding
#   threshold was not reached within the <= 5% audit path.
#
# Figure 2:
#   Contrasting influential-set persistence structures in a 1 x 2 panel.
#
#   Panel A:
#     Nested / stable-core influential-set evolution.
#
#   Panel B:
#     Non-nested / turnover influential-set evolution.
#
#   Visual encoding:
#     - blue tile       : observation belongs to the selected MIS at k
#     - orange diamond  : first entry into the MIS
#     - purple circle   : re-entry after previously leaving the MIS
#     - purple x        : exact non-nested MIS event at that deletion size
#     - dashed line     : first zero crossing of the audited coefficient
#
#   Observations are represented by their rank in order of first entry rather
#   than by arbitrary raw observation IDs.
#
#   The MIS membership colour is #2488BA.
#
# Figure 3:
#   NOT generated here.
#   The script exports the exact MIS observation IDs at the first zero crossing.
#   These IDs can later be joined to community/town/geographic metadata for a
#   genuine spatial-concentration figure.
#
# Visual design:
#   - soft / bright qualitative colours
#   - no plot titles; captions belong in LaTeX
#   - no red-green-only contrast
#   - PDF is the publication format
# ==============================================================================


# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "scales",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running scripts/90_output_application.R."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})


# ==============================================================================
# 1. Project paths
# ==============================================================================

resolve_project_root <- function() {
  
  candidates <- unique(c(
    normalizePath(
      ".",
      winslash = "/",
      mustWork = FALSE
    ),
    normalizePath(
      "..",
      winslash = "/",
      mustWork = FALSE
    )
  ))
  
  valid <- candidates[
    dir.exists(file.path(candidates, "scripts")) &
      dir.exists(file.path(candidates, "application"))
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

application_root <- file.path(
  project_root,
  "application"
)

output_root <- file.path(
  application_root,
  "00_cross_study"
)

fig_main_dir <- file.path(
  output_root,
  "figures",
  "main"
)

tab_main_dir <- file.path(
  output_root,
  "tables",
  "main"
)

data_dir <- file.path(
  output_root,
  "data"
)

diag_dir <- file.path(
  output_root,
  "diagnostics"
)

for (path in c(
  fig_main_dir,
  tab_main_dir,
  data_dir,
  diag_dir
)) {
  dir.create(
    path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==============================================================================
# 2. Study metadata
# ==============================================================================

study_meta <- data.frame(
  study_id = c(
    "01_RSPitHUMM",
    "02_LMDIVUiSL",
    "03_WbHS",
    "04_SBN-HIWTO",
    "05_TPCvIKT",
    "06_UHGtBtCEoaUWRP",
    "07_SwSIiUS",
    "08_IEoAtF",
    "09_StMEN",
    "10_TPoI"
  ),
  
  study_no = sprintf(
    "%02d",
    1:10
  ),
  
  study_label = c(
    "01  Mobile money",
    "02  Vaccine delivery",
    "03  Self-esteem",
    "04  HIV information",
    "05  Cash vs in-kind",
    "06  Workforce grant",
    "07  Sanitation spillovers",
    "08  Access to finance",
    "09  Immunization nudge",
    "10  Electricity information"
  ),
  
  stringsAsFactors = FALSE
)


study_ids <- study_meta$study_id


# ==============================================================================
# 3. Visual design
# ==============================================================================

# Soft but relatively vivid colours for the three threshold concepts.
COL_ATTENUATION <- "#6EA8FE"
COL_AMPLIFICATION <- "#F2B66D"
COL_CROSS_ZERO <- "#A78BFA"

threshold_colors <- c(
  "50% attenuation" = COL_ATTENUATION,
  "50% amplification" = COL_AMPLIFICATION,
  "Cross zero" = COL_CROSS_ZERO
)

threshold_shapes <- c(
  "50% attenuation" = 21,
  "50% amplification" = 22,
  "Cross zero" = 24
)


# Ten soft / bright study colours for Figure 2.

COL_GREY_DARK <- "#555555"
COL_GREY <- "#A7A7A7"
COL_GREY_LIGHT <- "#E2E2E2"
COL_GRID <- "#ECECEC"

SAVE_PNG_PREVIEWS <- FALSE
PNG_DPI <- 320

# ------------------------------------------------------------------------------
# Figure 2: contrasting influential-set persistence
# ------------------------------------------------------------------------------

# Panel A:
# A deliberately regular / nested example.
#
# Paper 08 produces a clean nested accumulation over the first 40 deletion
# sizes on the sign-reversing search path.
NESTED_STUDY <- "08_IEoAtF"
NESTED_PANEL_MAX_K <- 40L


# Panel B:
# Paper 04, Decrease direction, provides a readable non-nested example.
# Non-nested MIS events begin around k = 71 and the same path crosses zero
# at k = 95.

NON_NESTED_STUDY <- "04_SBN-HIWTO"
NON_NESTED_DIRECTION <- "Decrease"

NON_NESTED_PANEL_MIN_K <- 60L
NON_NESTED_PANEL_MAX_K <- 100L


# Exact requested MIS blue.
COL_PERSIST_SELECTED <- "#2488BA"

# Other soft / bright accents.
COL_PERSIST_ENTRY   <- "#F2B66D"
COL_PERSIST_REENTRY <- "#B79CED"
COL_PERSIST_EVENT   <- "#725AA8"
COL_PERSIST_CROSS   <- "#666666"

# ==============================================================================
# 4. Helpers
# ==============================================================================

read_required_csv <- function(
    study_id,
    filename
) {
  
  path <- file.path(
    application_root,
    study_id,
    filename
  )
  
  if (!file.exists(path)) {
    
    stop(
      "Missing required file:\n  ",
      path,
      "\n\nThis script expects the audit CSV files to exist under ",
      "application/<study>/."
    )
  }
  
  out <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  out
}


require_columns <- function(
    data,
    columns,
    source_name
) {
  
  missing <- setdiff(
    columns,
    names(data)
  )
  
  if (length(missing) > 0L) {
    stop(
      source_name,
      " is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
  
  invisible(TRUE)
}


safe_logical <- function(x) {
  
  if (is.logical(x)) {
    return(x)
  }
  
  x_chr <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  x_chr %in% c(
    "TRUE",
    "T",
    "1"
  )
}

as_nested_flag <- function(x) {
  
  if (is.logical(x)) {
    return(x)
  }
  
  x_chr <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  out <- rep(
    NA,
    length(x_chr)
  )
  
  out[
    x_chr %in% c(
      "TRUE",
      "T",
      "1"
    )
  ] <- TRUE
  
  out[
    x_chr %in% c(
      "FALSE",
      "F",
      "0"
    )
  ] <- FALSE
  
  out
}

first_threshold <- function(
    data,
    predicate,
    tie = c(
      "lower_ratio",
      "higher_ratio"
    )
) {
  
  tie <- match.arg(tie)
  
  keep <- predicate(
    data$normalized_ratio
  )
  
  keep[
    is.na(keep)
  ] <- FALSE
  
  hit <- data[
    keep,
    ,
    drop = FALSE
  ]
  
  if (nrow(hit) == 0L) {
    
    return(
      list(
        reached = FALSE,
        k = NA_integer_,
        fraction = NA_real_,
        direction = NA_character_,
        ratio = NA_real_,
        mis_ids = NA_character_
      )
    )
  }
  
  if (tie == "lower_ratio") {
    
    ord <- order(
      hit$removal_fraction,
      hit$k,
      hit$normalized_ratio
    )
    
  } else {
    
    ord <- order(
      hit$removal_fraction,
      hit$k,
      -hit$normalized_ratio
    )
  }
  
  row <- hit[
    ord[[1L]],
    ,
    drop = FALSE
  ]
  
  list(
    reached = TRUE,
    k = as.integer(row$k[[1L]]),
    fraction = as.numeric(row$removal_fraction[[1L]]),
    direction = as.character(row$direction[[1L]]),
    ratio = as.numeric(row$normalized_ratio[[1L]]),
    mis_ids = as.character(row$mis_ids[[1L]])
  )
}


format_beta <- function(x) {
  
  if (is.na(x)) {
    return("")
  }
  
  ax <- abs(x)
  
  if (ax >= 10) {
    return(
      sprintf(
        "%.2f",
        x
      )
    )
  }
  
  if (ax >= 1) {
    return(
      sprintf(
        "%.3f",
        x
      )
    )
  }
  
  if (ax >= 0.01) {
    return(
      sprintf(
        "%.4f",
        x
      )
    )
  }
  
  sprintf(
    "%.5f",
    x
  )
}


format_fraction_percent <- function(x) {
  
  if (is.na(x)) {
    return("Not reached")
  }
  
  pct <- 100 * x
  
  if (pct < 0.1) {
    
    return(
      paste0(
        sprintf(
          "%.3f",
          pct
        ),
        "\\%"
      )
    )
  }
  
  paste0(
    sprintf(
      "%.2f",
      pct
    ),
    "\\%"
  )
}


format_threshold <- function(
    fraction,
    k
) {
  
  if (
    is.na(fraction) ||
    is.na(k)
  ) {
    return("Not reached")
  }
  
  paste0(
    format_fraction_percent(
      fraction
    ),
    " (",
    as.integer(k),
    ")"
  )
}


save_plot <- function(
    plot,
    filename,
    width,
    height
) {
  
  pdf_file <- file.path(
    fig_main_dir,
    paste0(
      filename,
      ".pdf"
    )
  )
  
  ggsave(
    filename = pdf_file,
    plot = plot,
    width = width,
    height = height,
    units = "in"
  )
  
  if (SAVE_PNG_PREVIEWS) {
    
    png_file <- file.path(
      fig_main_dir,
      paste0(
        filename,
        ".png"
      )
    )
    
    ggsave(
      filename = png_file,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = PNG_DPI
    )
  }
  
  invisible(pdf_file)
}


# ==============================================================================
# 5. Read all ten audit outputs
# ==============================================================================

message(
  "Reading application audit outputs..."
)

path_list <- lapply(
  study_ids,
  function(id) {
    
    x <- read_required_csv(
      id,
      "audit_path.csv"
    )
    
    x$source_study <- id
    
    x
  }
)

baseline_list <- lapply(
  study_ids,
  function(id) {
    
    x <- read_required_csv(
      id,
      "audit_baseline.csv"
    )
    
    x$source_study <- id
    
    x
  }
)

summary_list <- lapply(
  study_ids,
  function(id) {
    
    x <- read_required_csv(
      id,
      "audit_summary.csv"
    )
    
    x$source_study <- id
    
    x
  }
)


audit_path <- bind_rows(
  path_list
)

audit_baseline <- bind_rows(
  baseline_list
)

audit_summary <- bind_rows(
  summary_list
)


# ==============================================================================
# 6. Validate input structure
# ==============================================================================

require_columns(
  audit_path,
  c(
    "study_id",
    "k",
    "removal_fraction",
    "direction",
    "beta_original",
    "beta_after",
    "mis_ids",
    "nested",
    "valid_refit"
  ),
  "Combined audit_path.csv"
)

require_columns(
  audit_baseline,
  c(
    "study_id",
    "n",
    "p",
    "beta_original"
  ),
  "Combined audit_baseline.csv"
)

require_columns(
  audit_summary,
  c(
    "study_id",
    "min_fraction_cross_zero"
  ),
  "Combined audit_summary.csv"
)


found_studies <- sort(
  unique(
    audit_path$study_id
  )
)

missing_studies <- setdiff(
  study_ids,
  found_studies
)

if (length(missing_studies) > 0L) {
  stop(
    "audit_path.csv files are missing study/studies: ",
    paste(
      missing_studies,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 7. Construct normalized paths
# ==============================================================================

audit_path$valid_refit_flag <- safe_logical(
  audit_path$valid_refit
)

n_before <- nrow(
  audit_path
)

audit_path_clean <- audit_path %>%
  filter(
    valid_refit_flag,
    is.finite(beta_original),
    is.finite(beta_after),
    is.finite(removal_fraction),
    beta_original != 0
  )


n_after <- nrow(
  audit_path_clean
)

if (n_after < n_before) {
  
  message(
    "Removed ",
    n_before - n_after,
    " invalid or unusable path rows."
  )
}


audit_path_clean <- audit_path_clean %>%
  mutate(
    normalized_ratio =
      beta_after / beta_original,
    
    nested_flag =
      as_nested_flag(
        nested
      ),
    
    path_type = case_when(
      
      beta_original > 0 &
        direction == "Increase" ~
        "Amplification search",
      
      beta_original > 0 &
        direction == "Decrease" ~
        "Attenuation / reversal search",
      
      beta_original < 0 &
        direction == "Decrease" ~
        "Amplification search",
      
      beta_original < 0 &
        direction == "Increase" ~
        "Attenuation / reversal search",
      
      TRUE ~
        as.character(direction)
    )
  ) %>%
  
  left_join(
    study_meta,
    by = "study_id"
  )


if (
  any(
    is.na(
      audit_path_clean$study_no
    )
  )
) {
  stop(
    "At least one study_id in audit_path.csv is not represented in study_meta."
  )
}


# Save the full normalized dataset used by the figures.
write.csv(
  audit_path_clean,
  file.path(
    data_dir,
    "05_cross_study_paths.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. Derive the three cross-study thresholds
# ==============================================================================

threshold_rows <- lapply(
  study_ids,
  function(id) {
    
    dat <- audit_path_clean %>%
      filter(
        study_id == id
      )
    
    if (nrow(dat) == 0L) {
      stop(
        "No valid path rows remain for ",
        id
      )
    }
    
    attenuation <- first_threshold(
      dat,
      predicate = function(r) {
        r <= 0.5
      },
      tie = "lower_ratio"
    )
    
    amplification <- first_threshold(
      dat,
      predicate = function(r) {
        r >= 1.5
      },
      tie = "higher_ratio"
    )
    
    cross_zero <- first_threshold(
      dat,
      predicate = function(r) {
        r <= 0
      },
      tie = "lower_ratio"
    )
    
    data.frame(
      study_id = id,
      
      max_audited_fraction =
        max(
          dat$removal_fraction,
          na.rm = TRUE
        ),
      
      attenuation_reached =
        attenuation$reached,
      
      attenuation_k =
        attenuation$k,
      
      attenuation_fraction =
        attenuation$fraction,
      
      attenuation_direction =
        attenuation$direction,
      
      attenuation_ratio =
        attenuation$ratio,
      
      amplification_reached =
        amplification$reached,
      
      amplification_k =
        amplification$k,
      
      amplification_fraction =
        amplification$fraction,
      
      amplification_direction =
        amplification$direction,
      
      amplification_ratio =
        amplification$ratio,
      
      cross_zero_reached =
        cross_zero$reached,
      
      cross_zero_k =
        cross_zero$k,
      
      cross_zero_fraction =
        cross_zero$fraction,
      
      cross_zero_direction =
        cross_zero$direction,
      
      cross_zero_ratio =
        cross_zero$ratio,
      
      cross_zero_mis_ids =
        cross_zero$mis_ids,
      
      stringsAsFactors = FALSE
    )
  }
)


thresholds <- bind_rows(
  threshold_rows
) %>%
  
  left_join(
    study_meta,
    by = "study_id"
  )


# ==============================================================================
# 9. Add baseline N, p and beta
# ==============================================================================

baseline_compact <- audit_baseline %>%
  group_by(
    study_id
  ) %>%
  summarise(
    N = first(n),
    p = first(p),
    beta_original_baseline =
      first(beta_original),
    .groups = "drop"
  )


thresholds <- thresholds %>%
  left_join(
    baseline_compact,
    by = "study_id"
  )


# Validate path beta against baseline beta.
path_beta_check <- audit_path_clean %>%
  group_by(
    study_id
  ) %>%
  summarise(
    beta_original_path =
      first(beta_original),
    .groups = "drop"
  )


thresholds <- thresholds %>%
  left_join(
    path_beta_check,
    by = "study_id"
  ) %>%
  
  mutate(
    beta_difference =
      beta_original_path -
      beta_original_baseline
  )


if (
  any(
    abs(
      thresholds$beta_difference
    ) > 1e-6,
    na.rm = TRUE
  )
) {
  
  warning(
    "At least one path beta_original differs from audit_baseline.csv ",
    "by more than 1e-6. Inspect 05_cross_study_thresholds.csv."
  )
}


# ==============================================================================
# 10. Validate zero-crossing calculation against existing audit_summary.csv
# ==============================================================================

summary_cross_check <- audit_summary %>%
  select(
    study_id,
    summary_cross_fraction =
      min_fraction_cross_zero
  )


thresholds <- thresholds %>%
  left_join(
    summary_cross_check,
    by = "study_id"
  ) %>%
  
  mutate(
    cross_fraction_difference =
      cross_zero_fraction -
      summary_cross_fraction
  )


bad_cross_check <- thresholds %>%
  filter(
    !is.na(cross_zero_fraction),
    !is.na(summary_cross_fraction),
    abs(cross_fraction_difference) > 1e-10
  )


if (nrow(bad_cross_check) > 0L) {
  
  warning(
    "The newly calculated zero-crossing threshold differs from ",
    "audit_summary.csv for at least one study. ",
    "Inspect 05_cross_study_thresholds.csv."
  )
}


write.csv(
  thresholds,
  file.path(
    data_dir,
    "05_cross_study_thresholds.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. Export zero-crossing MIS IDs for the future Figure 3
# ==============================================================================

cross_zero_ids <- thresholds %>%
  select(
    study_id,
    study_no,
    study_label,
    cross_zero_k,
    cross_zero_fraction,
    cross_zero_direction,
    cross_zero_ratio,
    cross_zero_mis_ids
  )


write.csv(
  cross_zero_ids,
  file.path(
    data_dir,
    "05_cross_zero_mis_ids.csv"
  ),
  row.names = FALSE
)


expand_id_row <- function(
    study_id,
    study_no,
    study_label,
    k,
    fraction,
    direction,
    ids
) {
  
  if (
    is.na(ids) ||
    !nzchar(
      trimws(ids)
    )
  ) {
    
    return(
      data.frame(
        study_id = character(),
        study_no = character(),
        study_label = character(),
        cross_zero_k = integer(),
        cross_zero_fraction = numeric(),
        cross_zero_direction = character(),
        observation_id = character(),
        stringsAsFactors = FALSE
      )
    )
  }
  
  id_vector <- strsplit(
    ids,
    split = ";",
    fixed = TRUE
  )[[1L]]
  
  data.frame(
    study_id = rep(
      study_id,
      length(id_vector)
    ),
    
    study_no = rep(
      study_no,
      length(id_vector)
    ),
    
    study_label = rep(
      study_label,
      length(id_vector)
    ),
    
    cross_zero_k = rep(
      k,
      length(id_vector)
    ),
    
    cross_zero_fraction = rep(
      fraction,
      length(id_vector)
    ),
    
    cross_zero_direction = rep(
      direction,
      length(id_vector)
    ),
    
    observation_id = trimws(
      id_vector
    ),
    
    stringsAsFactors = FALSE
  )
}


cross_zero_ids_long <- bind_rows(
  lapply(
    seq_len(
      nrow(
        cross_zero_ids
      )
    ),
    function(i) {
      
      row <- cross_zero_ids[
        i,
        ,
        drop = FALSE
      ]
      
      expand_id_row(
        study_id =
          row$study_id[[1L]],
        
        study_no =
          row$study_no[[1L]],
        
        study_label =
          row$study_label[[1L]],
        
        k =
          row$cross_zero_k[[1L]],
        
        fraction =
          row$cross_zero_fraction[[1L]],
        
        direction =
          row$cross_zero_direction[[1L]],
        
        ids =
          row$cross_zero_mis_ids[[1L]]
      )
    }
  )
)


write.csv(
  cross_zero_ids_long,
  file.path(
    data_dir,
    "05_cross_zero_mis_ids_long.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. Figure 1 data:
#     three sensitivity thresholds per study
# ==============================================================================

threshold_long <- bind_rows(
  
  thresholds %>%
    transmute(
      study_id,
      study_no,
      study_label,
      metric =
        "50% attenuation",
      reached =
        attenuation_reached,
      k =
        attenuation_k,
      fraction =
        attenuation_fraction
    ),
  
  thresholds %>%
    transmute(
      study_id,
      study_no,
      study_label,
      metric =
        "50% amplification",
      reached =
        amplification_reached,
      k =
        amplification_k,
      fraction =
        amplification_fraction
    ),
  
  thresholds %>%
    transmute(
      study_id,
      study_no,
      study_label,
      metric =
        "Cross zero",
      reached =
        cross_zero_reached,
      k =
        cross_zero_k,
      fraction =
        cross_zero_fraction
    )
)


threshold_long$metric <- factor(
  threshold_long$metric,
  levels = c(
    "50% attenuation",
    "50% amplification",
    "Cross zero"
  )
)


# Sort Figure 1 by the zero-crossing threshold.
# Studies without a zero crossing are placed at the end.
study_order <- thresholds %>%
  mutate(
    order_key = if_else(
      !is.na(
        cross_zero_fraction
      ),
      cross_zero_fraction,
      0.06 +
        if_else(
          !is.na(
            attenuation_fraction
          ),
          attenuation_fraction,
          0.06
        )
    )
  ) %>%
  arrange(
    order_key
  ) %>%
  pull(
    study_id
  )


plot_positions <- data.frame(
  study_id = study_order,
  
  study_y = rev(
    seq_along(
      study_order
    )
  ),
  
  stringsAsFactors = FALSE
)


threshold_long <- threshold_long %>%
  left_join(
    plot_positions,
    by = "study_id"
  )


metric_offsets <- c(
  "50% attenuation" = 0.22,
  "50% amplification" = 0,
  "Cross zero" = -0.22
)


# All audits run to approximately 5% of N.
# Place unachieved thresholds just to the right of the 5% line.
NOT_REACHED_X <- 0.052


threshold_long <- threshold_long %>%
  mutate(
    metric_offset =
      metric_offsets[
        as.character(metric)
      ],
    
    plot_y =
      study_y +
      metric_offset,
    
    plot_fraction =
      if_else(
        reached,
        fraction,
        NOT_REACHED_X
      )
  )


figure1_ranges <- threshold_long %>%
  group_by(
    study_id,
    study_y
  ) %>%
  summarise(
    x_min =
      min(
        plot_fraction,
        na.rm = TRUE
      ),
    x_max =
      max(
        plot_fraction,
        na.rm = TRUE
      ),
    .groups = "drop"
  )


figure1_labels <- study_meta %>%
  inner_join(
    plot_positions,
    by = "study_id"
  ) %>%
  arrange(
    study_y
  )


# ==============================================================================
# 13. Figure 1:
#     cross-study sensitivity thresholds
# ==============================================================================

fig1 <- ggplot() +
  
  geom_vline(
    xintercept = 0.05,
    colour = COL_GREY,
    linewidth = 0.45,
    linetype = "dashed"
  ) +
  
  geom_segment(
    data = figure1_ranges,
    aes(
      x = x_min,
      xend = x_max,
      y = study_y,
      yend = study_y
    ),
    colour = COL_GREY_LIGHT,
    linewidth = 0.8
  ) +
  
  geom_point(
    data = threshold_long %>%
      filter(
        reached
      ),
    
    aes(
      x = plot_fraction,
      y = plot_y,
      colour = metric,
      fill = metric,
      shape = metric
    ),
    
    size = 3.2,
    stroke = 0.75
  ) +
  
  geom_point(
    data = threshold_long %>%
      filter(
        !reached
      ),
    
    aes(
      x = plot_fraction,
      y = plot_y,
      colour = metric,
      shape = metric
    ),
    
    fill = "white",
    size = 3.2,
    stroke = 1
  ) +
  
  scale_colour_manual(
    values = threshold_colors,
    name = NULL
  ) +
  
  scale_fill_manual(
    values = threshold_colors,
    guide = "none"
  ) +
  
  scale_shape_manual(
    values = threshold_shapes,
    name = NULL
  ) +
  
  scale_x_continuous(
    breaks = seq(
      0,
      0.05,
      by = 0.01
    ),
    
    labels = label_percent(
      accuracy = 1
    ),
    
    limits = c(
      0,
      0.053
    ),
    
    expand = expansion(
      mult = c(
        0.01,
        0.01
      )
    )
  ) +
  
  scale_y_continuous(
    breaks =
      figure1_labels$study_y,
    
    labels =
      figure1_labels$study_label,
    
    expand = expansion(
      add = c(
        0.6,
        0.6
      )
    )
  ) +
  
  labs(
    x =
      "Fraction of estimation sample removed",
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 10.5
  ) +
  
  theme(
    panel.grid.minor =
      element_blank(),
    
    panel.grid.major.y =
      element_blank(),
    
    panel.grid.major.x =
      element_line(
        colour = COL_GRID,
        linewidth = 0.4
      ),
    
    axis.text.y =
      element_text(
        colour = COL_GREY_DARK,
        size = 9.3
      ),
    
    axis.text.x =
      element_text(
        colour = COL_GREY_DARK
      ),
    
    axis.title.x =
      element_text(
        margin = margin(
          t = 8
        )
      ),
    
    legend.position =
      "bottom",
    
    legend.direction =
      "horizontal",
    
    legend.box =
      "horizontal",
    
    plot.margin =
      margin(
        7,
        12,
        6,
        7
      )
  )


save_plot(
  fig1,
  filename =
    "05_fig1_sensitivity_thresholds",
  width = 7.4,
  height = 5.3
)


# Save exact Figure 1 plotting data.
write.csv(
  threshold_long,
  file.path(
    data_dir,
    "05_fig1_sensitivity_thresholds_data.csv"
  ),
  row.names = FALSE
)

# ==============================================================================
# 14. Figure 2:
#     nested versus non-nested influential-set persistence
# ==============================================================================

# ------------------------------------------------------------------------------
# 14.3 Build one persistence panel
# ------------------------------------------------------------------------------

build_persistence_panel <- function(
    study_id,
    direction,
    min_k,
    max_k,
    panel_role = c(
      "nested",
      "non_nested"
    ),
    show_y_axis = TRUE
) {
  
  panel_role <- match.arg(panel_role)
  
  
  # --------------------------------------------------------------------------
  # Full audit path for the selected study and direction
  # --------------------------------------------------------------------------
  
  available <- audit_path_clean %>%
    filter(
      .data$study_id == .env$study_id,
      .data$direction == .env$direction
    ) %>%
    arrange(
      k
    )
  
  
  if (nrow(available) == 0L) {
    stop(
      "No audit-path rows found for ",
      study_id,
      " in direction ",
      direction,
      "."
    )
  }
  
  
  max_available_k <- max(
    available$k,
    na.rm = TRUE
  )
  
  
  min_k <- max(
    1L,
    as.integer(min_k)
  )
  
  max_k <- min(
    as.integer(max_k),
    max_available_k
  )
  
  
  if (min_k > max_k) {
    stop(
      "Invalid persistence plotting window for ",
      study_id,
      ": min_k = ",
      min_k,
      ", max_k = ",
      max_k,
      "."
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Find zero crossing on THIS selected direction
  #
  # This is deliberately calculated from the path rather than from the
  # cross-study summary because the summary stores only the earliest crossing
  # across the two search directions.
  # --------------------------------------------------------------------------
  
  cross_rows <- available %>%
    filter(
      normalized_ratio <= 0
    ) %>%
    arrange(
      k
    )
  
  
  if (nrow(cross_rows) > 0L) {
    
    cross_zero_k <- as.integer(
      cross_rows$k[[1L]]
    )
    
    cross_zero_fraction <- as.numeric(
      cross_rows$removal_fraction[[1L]]
    )
    
  } else {
    
    cross_zero_k <- NA_integer_
    cross_zero_fraction <- NA_real_
  }
  
  
  # --------------------------------------------------------------------------
  # Use the complete history up to max_k to classify true first entry and
  # re-entry. Only afterwards restrict the displayed x-window.
  #
  # This prevents observations already present at min_k from being incorrectly
  # labelled as new entrants.
  # --------------------------------------------------------------------------
  
  history <- available %>%
    filter(
      k <= max_k
    )
  
  
  membership_history <- history %>%
    select(
      study_id,
      k,
      removal_fraction,
      direction,
      nested_flag,
      mis_ids
    ) %>%
    
    filter(
      !is.na(mis_ids),
      nzchar(
        trimws(
          mis_ids
        )
      )
    ) %>%
    
    tidyr::separate_rows(
      mis_ids,
      sep = ";"
    ) %>%
    
    mutate(
      observation_id =
        trimws(
          as.character(
            mis_ids
          )
        )
    ) %>%
    
    filter(
      nzchar(
        observation_id
      )
    ) %>%
    
    select(
      -mis_ids
    ) %>%
    
    distinct(
      study_id,
      k,
      removal_fraction,
      direction,
      nested_flag,
      observation_id
    )
  
  
  if (nrow(membership_history) == 0L) {
    stop(
      "No influential observation IDs available for ",
      study_id,
      "."
    )
  }
  
  
  # --------------------------------------------------------------------------
  # True membership history
  # --------------------------------------------------------------------------
  
  membership_history <- membership_history %>%
    arrange(
      observation_id,
      k
    ) %>%
    
    group_by(
      observation_id
    ) %>%
    
    mutate(
      previous_selected_k =
        lag(
          k
        ),
      
      membership_event =
        case_when(
          
          row_number() == 1L ~
            "First entry",
          
          k - previous_selected_k > 1L ~
            "Re-entry",
          
          TRUE ~
            "Retained"
        )
    ) %>%
    
    ungroup()
  
  
  # --------------------------------------------------------------------------
  # Restrict to the visible plotting window
  # --------------------------------------------------------------------------
  
  membership <- membership_history %>%
    filter(
      k >= min_k,
      k <= max_k
    )
  
  
  visible_ids <- membership %>%
    distinct(
      observation_id
    ) %>%
    pull(
      observation_id
    )
  
  
  observation_summary <- membership_history %>%
    filter(
      observation_id %in% visible_ids
    ) %>%
    
    group_by(
      observation_id
    ) %>%
    
    summarise(
      first_k =
        min(
          k
        ),
      
      last_k =
        max(
          k
        ),
      
      n_selected =
        n_distinct(
          k
        ),
      
      .groups = "drop"
    ) %>%
    
    arrange(
      first_k,
      desc(
        n_selected
      ),
      observation_id
    ) %>%
    
    mutate(
      observation_order =
        row_number()
    )
  
  
  membership <- membership %>%
    left_join(
      observation_summary %>%
        select(
          observation_id,
          observation_order,
          first_k
        ),
      by = "observation_id"
    )
  
  
  n_observations <- nrow(
    observation_summary
  )
  
  
  # --------------------------------------------------------------------------
  # Membership events visible in the selected window
  # --------------------------------------------------------------------------
  
  first_entry <- membership %>%
    filter(
      membership_event == "First entry"
    )
  
  
  reentry <- membership %>%
    filter(
      membership_event == "Re-entry"
    )
  
  
  non_nested_steps <- available %>%
    filter(
      k >= min_k,
      k <= max_k,
      !is.na(
        nested_flag
      ),
      nested_flag == FALSE
    ) %>%
    
    distinct(
      k
    ) %>%
    
    mutate(
      marker_y =
        n_observations +
        0.75
    )
  
  
  # --------------------------------------------------------------------------
  # Save panel data
  # --------------------------------------------------------------------------
  
  write.csv(
    membership,
    file.path(
      data_dir,
      paste0(
        "05_fig2_",
        panel_role,
        "_",
        study_id,
        "_membership.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # --------------------------------------------------------------------------
  # Plot
  # --------------------------------------------------------------------------
  
  p <- ggplot(
    membership,
    aes(
      x = k,
      y = observation_order
    )
  ) +
    
    geom_tile(
      fill = COL_PERSIST_SELECTED,
      width = 0.90,
      height = 0.82,
      alpha = 0.72
    ) +
    
    geom_point(
      data = first_entry,
      aes(
        x = k,
        y = observation_order
      ),
      inherit.aes = FALSE,
      shape = 23,
      size = 1.85,
      stroke = 0.45,
      fill = COL_PERSIST_ENTRY,
      colour = "white"
    ) +
    
    geom_point(
      data = reentry,
      aes(
        x = k,
        y = observation_order
      ),
      inherit.aes = FALSE,
      shape = 21,
      size = 1.75,
      stroke = 0.45,
      fill = COL_PERSIST_REENTRY,
      colour = "white"
    )
  
  
  # Zero-crossing line is shown only when it lies inside the displayed window.
  if (
    !is.na(cross_zero_k) &&
    cross_zero_k >= min_k &&
    cross_zero_k <= max_k
  ) {
    
    p <- p +
      
      geom_vline(
        xintercept = cross_zero_k,
        colour = COL_PERSIST_CROSS,
        linewidth = 0.55,
        linetype = "dashed"
      )
  }
  
  
  # Exact non-nested MIS events.
  if (nrow(non_nested_steps) > 0L) {
    
    p <- p +
      
      geom_point(
        data = non_nested_steps,
        aes(
          x = k,
          y = marker_y
        ),
        inherit.aes = FALSE,
        shape = 4,
        size = 2.1,
        stroke = 0.8,
        colour = COL_PERSIST_EVENT
      )
  }
  
  
  y_breaks <- unique(
    round(
      pretty(
        c(
          1,
          n_observations
        ),
        n = 6
      )
    )
  )
  
  y_breaks <- y_breaks[
    y_breaks >= 1 &
      y_breaks <= n_observations
  ]
  
  
  p <- p +
    
    scale_x_continuous(
      breaks = pretty(
        c(
          min_k,
          max_k
        ),
        n = 6
      ),
      
      limits = c(
        min_k - 0.5,
        max_k + 0.5
      ),
      
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    
    scale_y_reverse(
      breaks = y_breaks,
      
      limits = c(
        n_observations +
          1.25,
        0.5
      ),
      
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    
    labs(
      subtitle =
        if (panel_role == "nested") {
          "Stable-core accumulation"
        } else {
          "Influential-set turnover"
        },
      
      x =
        "Deletion-set size k",
      
      y =
        if (show_y_axis) {
          "Observation rank by first entry"
        } else {
          NULL
        }
    ) +
    
    theme_minimal(
      base_size = 10.2
    ) +
    
    theme(
      panel.grid =
        element_blank(),
      
      plot.subtitle =
        element_text(
          face = "bold",
          size = 10.5,
          colour = COL_GREY_DARK,
          margin = margin(
            b = 7
          )
        ),
      
      axis.text =
        element_text(
          colour = COL_GREY_DARK
        ),
      
      axis.text.y =
        if (show_y_axis) {
          element_text(
            size = 8
          )
        } else {
          element_blank()
        },
      
      axis.ticks.y =
        if (show_y_axis) {
          element_line(
            colour = COL_GREY_LIGHT,
            linewidth = 0.3
          )
        } else {
          element_blank()
        },
      
      axis.title.x =
        element_text(
          margin = margin(
            t = 7
          )
        ),
      
      axis.title.y =
        element_text(
          margin = margin(
            r = 7
          )
        ),
      
      plot.margin =
        margin(
          6,
          5,
          5,
          5
        )
    )
  
  
  list(
    plot = p,
    
    study_id =
      study_id,
    
    direction =
      direction,
    
    min_k =
      min_k,
    
    max_k =
      max_k,
    
    cross_zero_k =
      cross_zero_k,
    
    cross_zero_fraction =
      cross_zero_fraction,
    
    n_observations =
      n_observations,
    
    n_non_nested =
      nrow(
        non_nested_steps
      ),
    
    n_reentries =
      nrow(
        reentry
      )
  )
}

# ------------------------------------------------------------------------------
# 14.4 Build Panel A: nested / stable-core path
# ------------------------------------------------------------------------------

nested_panel <- build_persistence_panel(
  study_id =
    NESTED_STUDY,
  
  direction =
    "Increase",
  
  min_k =
    1L,
  
  max_k =
    NESTED_PANEL_MAX_K,
  
  panel_role =
    "nested",
  
  show_y_axis =
    TRUE
)


# ------------------------------------------------------------------------------
# 14.5 Build Panel B: non-nested / turnover path
# ------------------------------------------------------------------------------

non_nested_panel <- build_persistence_panel(
  study_id =
    NON_NESTED_STUDY,
  
  direction =
    NON_NESTED_DIRECTION,
  
  min_k =
    NON_NESTED_PANEL_MIN_K,
  
  max_k =
    NON_NESTED_PANEL_MAX_K,
  
  panel_role =
    "non_nested",
  
  show_y_axis =
    FALSE
)


# ------------------------------------------------------------------------------
# Shared visual key for Figure 2
# ------------------------------------------------------------------------------

legend_plot <- ggplot() +
  
  # Selected MIS membership.
  annotate(
    "rect",
    xmin = 0.55,
    xmax = 0.85,
    ymin = 0.36,
    ymax = 0.64,
    fill = COL_PERSIST_SELECTED,
    alpha = 0.72
  ) +
  
  annotate(
    "text",
    x = 0.95,
    y = 0.50,
    label = "Selected in MIS",
    hjust = 0,
    size = 3.15,
    colour = COL_GREY_DARK
  ) +
  
  # First entry.
  annotate(
    "point",
    x = 2.45,
    y = 0.50,
    shape = 23,
    size = 3.1,
    stroke = 0.55,
    fill = COL_PERSIST_ENTRY,
    colour = "white"
  ) +
  
  annotate(
    "text",
    x = 2.62,
    y = 0.50,
    label = "First entry",
    hjust = 0,
    size = 3.15,
    colour = COL_GREY_DARK
  ) +
  
  # Re-entry.
  annotate(
    "point",
    x = 3.75,
    y = 0.50,
    shape = 21,
    size = 3.0,
    stroke = 0.55,
    fill = COL_PERSIST_REENTRY,
    colour = "white"
  ) +
  
  annotate(
    "text",
    x = 3.92,
    y = 0.50,
    label = "Re-entry",
    hjust = 0,
    size = 3.15,
    colour = COL_GREY_DARK
  ) +
  
  # Non-nested step.
  annotate(
    "point",
    x = 4.90,
    y = 0.50,
    shape = 4,
    size = 3.0,
    stroke = 0.85,
    colour = COL_PERSIST_EVENT
  ) +
  
  annotate(
    "text",
    x = 5.08,
    y = 0.50,
    label = "Non-nested step",
    hjust = 0,
    size = 3.15,
    colour = COL_GREY_DARK
  ) +
  
  # Zero crossing.
  annotate(
    "segment",
    x = 6.75,
    xend = 6.75,
    y = 0.32,
    yend = 0.68,
    colour = COL_PERSIST_CROSS,
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  
  annotate(
    "text",
    x = 6.90,
    y = 0.50,
    label = "Zero crossing",
    hjust = 0,
    size = 3.15,
    colour = COL_GREY_DARK
  ) +
  
  coord_cartesian(
    xlim = c(
      0.4,
      8.0
    ),
    ylim = c(
      0.25,
      0.75
    ),
    clip = "off"
  ) +
  
  theme_void() +
  
  theme(
    plot.margin =
      margin(
        0,
        4,
        0,
        4
      )
  )

# ------------------------------------------------------------------------------
# 14.6 Combine into a 1 x 2 publication figure
# ------------------------------------------------------------------------------

fig2_panels <- (
  nested_panel$plot |
    non_nested_panel$plot
) +
  
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  
  theme(
    plot.tag =
      element_text(
        face = "bold",
        size = 12
      )
  )


fig2 <- (
  fig2_panels /
    legend_plot
) +
  
  patchwork::plot_layout(
    heights = c(
      1,
      0.09
    )
  )


save_plot(
  fig2,
  filename =
    "05_fig2_influential_set_persistence",
  width = 7.6,
  height = 5.0
)


# ------------------------------------------------------------------------------
# 14.7 Figure 2 diagnostics
# ------------------------------------------------------------------------------

fig2_diagnostics <- data.frame(
  panel = c(
    "A",
    "B"
  ),
  
  role = c(
    "Nested / stable-core",
    "Non-nested / turnover"
  ),
  
  study_id = c(
    nested_panel$study_id,
    non_nested_panel$study_id
  ),
  
  direction = c(
    nested_panel$direction,
    non_nested_panel$direction
  ),
  
  cross_zero_k = c(
    nested_panel$cross_zero_k,
    non_nested_panel$cross_zero_k
  ),
  
  cross_zero_fraction = c(
    nested_panel$cross_zero_fraction,
    non_nested_panel$cross_zero_fraction
  ),
  
  plotted_min_k = c(
    nested_panel$min_k,
    non_nested_panel$min_k
  ),
  
  plotted_max_k = c(
    nested_panel$max_k,
    non_nested_panel$max_k
  ),
  
  n_observations = c(
    nested_panel$n_observations,
    non_nested_panel$n_observations
  ),
  
  n_non_nested_events = c(
    nested_panel$n_non_nested,
    non_nested_panel$n_non_nested
  ),
  
  n_reentries = c(
    nested_panel$n_reentries,
    non_nested_panel$n_reentries
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  fig2_diagnostics,
  file.path(
    diag_dir,
    "05_fig2_persistence_panel_diagnostics.csv"
  ),
  row.names = FALSE
)


message("")
message("Figure 2 panel selection:")
message(
  "  Panel A: ",
  nested_panel$study_id,
  " | non-nested events shown = ",
  nested_panel$n_non_nested
)

message(
  "  Panel B: ",
  non_nested_panel$study_id,
  " | non-nested events shown = ",
  non_nested_panel$n_non_nested,
  " | re-entry events = ",
  non_nested_panel$n_reentries
)

# ==============================================================================
# 17. Main Table:
#     cross-study sensitivity thresholds
# ==============================================================================

table_data <- thresholds %>%
  arrange(
    study_no
  ) %>%
  
  transmute(
    study_id,
    study_no,
    Study =
      study_label,
    
    N =
      as.integer(N),
    
    p =
      as.integer(p),
    
    beta_original =
      beta_original_baseline,
    
    attenuation_k,
    attenuation_fraction,
    
    amplification_k,
    amplification_fraction,
    
    cross_zero_k,
    cross_zero_fraction
  )


write.csv(
  table_data,
  file.path(
    tab_main_dir,
    "05_tab1_cross_study_sensitivity.csv"
  ),
  row.names = FALSE
)


table_rows <- vapply(
  seq_len(
    nrow(
      table_data
    )
  ),
  function(i) {
    
    row <- table_data[
      i,
      ,
      drop = FALSE
    ]
    
    paste0(
      row$Study[[1L]],
      " & ",
      format(
        row$N[[1L]],
        big.mark = ",",
        scientific = FALSE
      ),
      " & ",
      format_beta(
        row$beta_original[[1L]]
      ),
      " & ",
      format_threshold(
        row$attenuation_fraction[[1L]],
        row$attenuation_k[[1L]]
      ),
      " & ",
      format_threshold(
        row$amplification_fraction[[1L]],
        row$amplification_k[[1L]]
      ),
      " & ",
      format_threshold(
        row$cross_zero_fraction[[1L]],
        row$cross_zero_k[[1L]]
      ),
      " \\\\"
    )
  },
  character(1)
)


table_tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  
  paste0(
    "\\caption{",
    "Cross-study target-coefficient sensitivity thresholds.",
    "}"
  ),
  
  "\\label{tab:application-sensitivity}",
  
  "\\small",
  
  "\\resizebox{\\textwidth}{!}{%",
  
  "\\begin{tabular}{lrrrrr}",
  
  "\\toprule",
  
  paste0(
    "Study",
    " & $N$",
    " & $\\widehat{\\beta}_{\\mathrm{full}}$",
    " & 50\\% attenuation",
    " & 50\\% amplification",
    " & Cross zero",
    " \\\\"
  ),
  
  "\\midrule",
  
  table_rows,
  
  "\\bottomrule",
  
  "\\end{tabular}%",
  
  "}",
  
  "\\begin{minipage}{0.98\\textwidth}",
  
  "\\footnotesize",
  
  paste0(
    "\\textit{Notes:} ",
    "Threshold entries report the fraction of the estimation sample removed, ",
    "with the corresponding number of deleted observations in parentheses. ",
    "The normalized coefficient is ",
    "$R_k=\\widehat{\\beta}_{-S_k}/\\widehat{\\beta}_{\\mathrm{full}}$. ",
    "A 50\\% attenuation is first reached when $R_k\\leq0.5$; ",
    "a 50\\% amplification is first reached when $R_k\\geq1.5$; ",
    "and the coefficient crosses zero when $R_k\\leq0$. ",
    "``Not reached'' indicates that the threshold is not attained within the ",
    "study's prespecified deletion budget of at most approximately 5\\% of ",
    "the estimation sample."
  ),
  
  "\\end{minipage}",
  
  "\\end{table}"
)


writeLines(
  table_tex,
  con = file.path(
    tab_main_dir,
    "05_tab1_cross_study_sensitivity.tex"
  )
)


# ==============================================================================
# 18. Figure 3 preparation note
# ==============================================================================

figure3_note <- c(
  
  "Figure 3 is intentionally not generated by scripts/90_output_application.R.",
  
  "",
  
  "Reason:",
  "The audit outputs identify the influential observations but do not by",
  "themselves provide geographic coordinates or interpretable location labels.",
  
  "",
  
  "Prepared inputs:",
  "  data/05_cross_zero_mis_ids.csv",
  "  data/05_cross_zero_mis_ids_long.csv",
  
  "",
  
  "Recommended next step:",
  "Choose the most informative application, join observation_id back to the",
  "paper-specific analysis sample, and retain meaningful location information",
  "such as community, town, village, district, latitude, or longitude.",
  
  "",
  
  "Do not construct a geographic map from arbitrary numeric cluster IDs alone.",
  
  "",
  
  "Potential Figure 3 quantity:",
  "MIS enrichment for geographic unit g =",
  "  (share of sign-reversing MIS observations in g) /",
  "  (share of the estimation sample in g)."
)


writeLines(
  figure3_note,
  con = file.path(
    diag_dir,
    "05_figure3_status.txt"
  )
)


# ==============================================================================
# 19. Diagnostics and manifest
# ==============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    diag_dir,
    "session_info.txt"
  )
)


manifest <- list.files(
  output_root,
  recursive = TRUE,
  full.names = FALSE
)


writeLines(
  manifest,
  con = file.path(
    diag_dir,
    "manifest.txt"
  )
)


# ==============================================================================
# 20. Console summary
# ==============================================================================

message("")
message("Section 5 cross-study outputs generated successfully.")
message("")
message("Output root:")
message("  ", output_root)
message("")
message("Main figures:")
message("  figures/main/05_fig1_sensitivity_thresholds.pdf")
message("  figures/main/05_fig2_influential_set_persistence.pdf")
message("")
message("Main table:")
message("  tables/main/05_tab1_cross_study_sensitivity.tex")
message("")
message("Figure 3 preparation:")
message("  data/05_cross_zero_mis_ids_long.csv")
message("  diagnostics/05_figure3_status.txt")
message("")