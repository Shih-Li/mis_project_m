# ==============================================================================
# File: R/sim_misspec_engine.R
# Purpose:
#   One Monte Carlo iteration for Script 05.
#
# Main optimisation:
#   - Generate each dataset ONCE
#   - Fit full OLS ONCE per model state
#   - Fit full MM ONCE per model state
#   - Prepare FWL/MIS problem ONCE per model state
#   - Compute DFBETAS/Cook/leverage ONCE per model state
#   - Reuse all of the above for every k/n
#
# Required sourced files:
#   R/diagnostics_classical.R
#   R/dgp_misspec_factory.R
#   R/mis_sensitivity.R
#   R/misspec_metrics.R
# ==============================================================================


# ==============================================================================
# 1. Locate target coefficient in model matrix
# ==============================================================================

get_target_pos05 <- function(
    model,
    target_var = "x"
) {
  
  mm <- stats::model.matrix(model)
  
  pos <- match(
    target_var,
    colnames(mm)
  )
  
  if (is.na(pos)) {
    stop(
      "Target variable not found in model matrix: ",
      target_var
    )
  }
  
  pos
}


# ==============================================================================
# 2. Extract target coefficient + SE from an ALREADY fitted lm object
#
# This avoids fitting the full OLS model twice.
# ==============================================================================

extract_target_lm05 <- function(
    model,
    target_var = "x"
) {
  
  out <- tryCatch({
    
    sm <- summary(model)$coefficients
    
    if (!target_var %in% rownames(sm)) {
      
      return(
        c(
          coef = NA_real_,
          se = NA_real_
        )
      )
    }
    
    c(
      coef = unname(
        sm[target_var, "Estimate"]
      ),
      se = unname(
        sm[target_var, "Std. Error"]
      )
    )
    
  }, error = function(e) {
    
    c(
      coef = NA_real_,
      se = NA_real_
    )
    
  })
  
  out
}


# ==============================================================================
# 3. Construct one result row
# ==============================================================================

make_method_row05 <- function(
    base_meta,
    model_state,
    diagnostic,
    estimator,
    fit,
    full_fit,
    abnormality,
    selected_idx,
    affected_idx,
    runtime = NA_real_
) {
  
  # --------------------------------------------------------------------------
  # Coefficient movement relative to full OLS
  # --------------------------------------------------------------------------
  
  delta <- if (
    is.finite(fit["coef"]) &&
    is.finite(full_fit["coef"])
  ) {
    
    abs(
      fit["coef"] -
        full_fit["coef"]
    )
    
  } else {
    
    NA_real_
    
  }
  
  
  # --------------------------------------------------------------------------
  # Standardise movement using SE of full OLS coefficient
  # --------------------------------------------------------------------------
  
  standardized_delta <- if (
    is.finite(delta) &&
    is.finite(full_fit["se"]) &&
    full_fit["se"] > 0
  ) {
    
    delta /
      full_fit["se"]
    
  } else {
    
    NA_real_
    
  }
  
  
  # --------------------------------------------------------------------------
  # Estimation accuracy
  # --------------------------------------------------------------------------
  
  em <- compute_estimation_metrics05(
    fit["coef"],
    fit["se"],
    base_meta$true_beta
  )
  
  # --------------------------------------------------------------------------
  # Additional detailed metrics
  # --------------------------------------------------------------------------
  
  coef_now <-
    as.numeric(
      fit["coef"]
    )[1L]
  
  
  full_coef <-
    as.numeric(
      full_fit["coef"]
    )[1L]
  
  
  full_se <-
    as.numeric(
      full_fit["se"]
    )[1L]
  
  
  # Signed coefficient movement
  signed_delta_beta <-
    if (
      is.finite(coef_now) &&
      is.finite(full_coef)
    ) {
      
      coef_now -
        full_coef
      
    } else {
      
      NA_real_
    }
  
  
  # Signed standardized movement
  signed_standardized_delta <-
    if (
      is.finite(signed_delta_beta) &&
      is.finite(full_se) &&
      full_se >
      .Machine$double.eps
    ) {
      
      signed_delta_beta /
        full_se
      
    } else {
      
      NA_real_
    }
  
  
  # Relative movement against magnitude of original OLS coefficient
  relative_delta_beta <-
    if (
      is.finite(signed_delta_beta) &&
      is.finite(full_coef) &&
      abs(full_coef) >
      sqrt(.Machine$double.eps)
    ) {
      
      signed_delta_beta /
        abs(full_coef)
      
    } else {
      
      NA_real_
    }
  
  
  # Unique selected observations
  selected_unique <-
    if (
      is.null(selected_idx)
    ) {
      
      integer(0)
      
    } else {
      
      unique(
        as.integer(
          selected_idx
        )
      )
    }
  
  
  # Genuine affected set exists only for subgroup-type DGPs
  affected_unique <-
    if (
      is.null(affected_idx) ||
      length(affected_idx) == 0L
    ) {
      
      NULL
      
    } else {
      
      unique(
        as.integer(
          affected_idx
        )
      )
    }
  
  
  affected_n <-
    if (
      is.null(affected_unique)
    ) {
      
      NA_integer_
      
    } else {
      
      length(
        affected_unique
      )
    }
  
  
  affected_fraction <-
    if (
      is.finite(affected_n) &&
      is.finite(base_meta$n) &&
      base_meta$n > 0
    ) {
      
      affected_n /
        base_meta$n
      
    } else {
      
      NA_real_
    }
  
  
  overlap_n <-
    if (
      is.null(affected_unique)
    ) {
      
      NA_integer_
      
    } else {
      
      length(
        intersect(
          selected_unique,
          affected_unique
        )
      )
    }
  
  
  # Selection enrichment relative to random selection.
  #
  # Example:
  # affected_fraction = 0.25
  # precision = 0.75
  # selection_lift = 3
  #
  # => selected observations are 3x more likely to belong to the
  #    genuine affected group than random selection.
  selection_lift <-
    if (
      !is.null(affected_unique) &&
      length(selected_unique) > 0L &&
      is.finite(affected_fraction) &&
      affected_fraction > 0
    ) {
      
      selected_precision <-
        overlap_n /
        length(selected_unique)
      
      selected_precision /
        affected_fraction
      
    } else {
      
      NA_real_
    }
  
  # --------------------------------------------------------------------------
  # Output
  # --------------------------------------------------------------------------
  
  data.frame(
    
    iter =
      base_meta$iter,
    
    n =
      base_meta$n,
    
    true_beta =
      base_meta$true_beta,
    
    draw_seed =
      base_meta$draw_seed,
    
    x_type =
      base_meta$x_type,
    
    error_type =
      base_meta$error_type,
    
    scenario =
      base_meta$scenario,
    
    severity_level =
      base_meta$severity_level,
    
    severity_target =
      base_meta$severity_target,
    
    severity_realized =
      base_meta$severity_realized,
    
    raw_misspec_coef =
      base_meta$raw_misspec_coef,
    
    k_fraction =
      base_meta$k_fraction,
    
    k =
      base_meta$k,
    
    model_state =
      model_state,
    
    diagnostic =
      diagnostic,
    
    estimator =
      estimator,
    
    coef =
      unname(
        fit["coef"]
      ),
    
    se =
      unname(
        fit["se"]
      ),
    
    signed_delta_beta =
      signed_delta_beta,
    
    signed_standardized_delta =
      signed_standardized_delta,
    
    relative_delta_beta =
      relative_delta_beta,
    
    delta_beta =
      unname(
        delta
      ),
    
    standardized_delta =
      unname(
        standardized_delta
      ),
    
    outlier_rate =
      abnormality$rate,
    
    n_abnormal =
      abnormality$n,
    
    residual_mad =
      abnormality$residual_mad,
    
    q95_abs_std_resid =
      abnormality$q95_abs_std_resid,
    
    max_abs_std_resid =
      abnormality$max_abs_std_resid,
    
    max_leverage =
      abnormality$max_leverage,
    
    condition_number =
      abnormality$condition_number,
    
    overlap_recall =
      compute_selection_overlap05(
        selected_idx,
        affected_idx
      ),
    
    overlap_precision =
      compute_precision05(
        selected_idx,
        affected_idx
      ),
    
    affected_n =
      affected_n,
    
    affected_fraction =
      affected_fraction,
    
    overlap_n =
      overlap_n,
    
    selection_lift =
      selection_lift,
    
    selected_size =
      length(selected_idx),
    
    bias =
      unname(
        em["bias"]
      ),
    
    abs_bias =
      unname(
        em["abs_bias"]
      ),
    
    coverage =
      unname(
        em["coverage"]
      ),
    
    runtime =
      runtime,
    
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 4. Process ONE model state using ALL k fractions
#
# Examples:
#
#   state_name = "wrong"
#   state_name = "correct"
#
# Full model quantities are calculated ONCE and reused for:
#
#   k/n = 1%
#   k/n = 2.5%
#   k/n = 5%
#   k/n = 10%
# ==============================================================================

run_state05 <- function(
    
  data,
  
  formula,
  
  state_name,
  
  k_fractions,
  
  affected_idx,
  
  base_meta,
  
  target_var = "x",
  
  use_cpp = TRUE
  
) {
  
  # ============================================================================
  # 4.1 Fit full OLS ONCE
  # ============================================================================
  
  t_ols <- proc.time()[3L]
  
  
  mod <- tryCatch(
    
    stats::lm(
      formula,
      data = data
    ),
    
    error = function(e) {
      
      stop(
        sprintf(
          paste0(
            "Full OLS failed: ",
            "scenario=%s, state=%s, n=%d. ",
            "Original error: %s"
          ),
          base_meta$scenario,
          state_name,
          base_meta$n,
          conditionMessage(e)
        )
      )
    }
  )
  
  
  full_ols <- extract_target_lm05(
    mod,
    target_var
  )
  
  
  full_ols_runtime <-
    proc.time()[3L] -
    t_ols
  
  
  # ============================================================================
  # 4.2 Fit full MM ONCE
  # ============================================================================
  
  t_mm_full <- proc.time()[3L]
  
  
  full_mm <- fit_target_mm05(
    formula,
    data,
    target_var
  )
  
  
  full_mm_runtime <-
    proc.time()[3L] -
    t_mm_full
  
  
  # ============================================================================
  # 4.3 Residual abnormality ONCE
  # ============================================================================
  
  abnormality <-
    compute_residual_abnormality05(
      mod
    )
  
  
  # ============================================================================
  # 4.4 Find target coefficient position ONCE
  # ============================================================================
  
  target_pos <-
    get_target_pos05(
      mod,
      target_var
    )
  
  
  # ============================================================================
  # 4.5 Prepare FWL / MIS quantities ONCE
  #
  # This avoids repeating:
  #
  #   model.matrix()
  #   QR decomposition
  #   FWL residualisation
  #
  # for every k and every MIS direction.
  # ============================================================================
  
  mis_problem <-
    prepare_mis_problem05(
      
      mod,
      
      pos =
        target_pos
    )
  
  
  # ============================================================================
  # 4.6 Compute classical diagnostics ONCE
  #
  # get_all_classical() returns:
  #
  #   id
  #   leverage
  #   cooks_d
  #   dfbetas_target
  # ============================================================================
  
  classical <-
    get_all_classical(
      mod,
      target_var
    )
  
  # --------------------------------------------------------------------------
  # Additional state-level diagnostics
  #
  # These are computed once per fitted model state and reused for all methods/k.
  # --------------------------------------------------------------------------
  
  abnormality$max_leverage <-
    if (
      length(classical$leverage) > 0L &&
      any(is.finite(classical$leverage))
    ) {
      
      max(
        classical$leverage,
        na.rm = TRUE
      )
      
    } else {
      
      NA_real_
    }
  
  
  abnormality$condition_number <-
    tryCatch(
      {
        
        X_model <-
          stats::model.matrix(
            mod
          )
        
        as.numeric(
          kappa(
            X_model,
            exact = FALSE
          )
        )
        
      },
      error = function(e) {
        NA_real_
      }
    )
  
  
  # --------------------------------------------------------------------------
  # Rank DFBETAS ONCE
  # --------------------------------------------------------------------------
  
  rank_dfb <-
    classical$id[
      order(
        abs(
          classical$dfbetas_target
        ),
        decreasing = TRUE
      )
    ]
  
  
  # --------------------------------------------------------------------------
  # Rank Cook's D ONCE
  # --------------------------------------------------------------------------
  
  rank_cook <-
    classical$id[
      order(
        classical$cooks_d,
        decreasing = TRUE
      )
    ]
  
  
  # --------------------------------------------------------------------------
  # Rank leverage ONCE
  # --------------------------------------------------------------------------
  
  rank_lev <-
    classical$id[
      order(
        classical$leverage,
        decreasing = TRUE
      )
    ]
  
  
  # ============================================================================
  # 4.7 Result storage
  # ============================================================================
  
  rows <- list()
  
  
  # ============================================================================
  # 4.8 Evaluate ALL k values using SAME dataset/model
  # ============================================================================
  
  for (
    k_fraction
    in
    k_fractions
  ) {
    
    # --------------------------------------------------------------------------
    # Convert fraction into integer k
    # --------------------------------------------------------------------------
    
    k <- max(
      
      1L,
      
      min(
        
        nrow(data) - 2L,
        
        floor(
          nrow(data) *
            k_fraction
        )
      )
    )
    
    
    # --------------------------------------------------------------------------
    # Add k-specific metadata
    # --------------------------------------------------------------------------
    
    meta <-
      base_meta
    
    meta$k_fraction <-
      k_fraction
    
    meta$k <-
      k
    
    
    # ==========================================================================
    # 4.9 Full OLS baseline
    #
    # IMPORTANT:
    # This does NOT refit OLS.
    #
    # We simply repeat the already calculated result as a row for this k so
    # later comparisons remain grouped by k_fraction.
    # ==========================================================================
    
    rows[[
      length(rows) + 1L
    ]] <-
      make_method_row05(
        
        base_meta =
          meta,
        
        model_state =
          state_name,
        
        diagnostic =
          "none",
        
        estimator =
          "OLS",
        
        fit =
          full_ols,
        
        full_fit =
          full_ols,
        
        abnormality =
          abnormality,
        
        selected_idx =
          integer(0),
        
        affected_idx =
          affected_idx,
        
        runtime =
          full_ols_runtime
      )
    
    
    # ==========================================================================
    # 4.10 Full MM baseline
    #
    # Also NOT refitted for each k.
    # ==========================================================================
    
    rows[[
      length(rows) + 1L
    ]] <-
      make_method_row05(
        
        base_meta =
          meta,
        
        model_state =
          state_name,
        
        diagnostic =
          "none",
        
        estimator =
          "MM",
        
        fit =
          full_mm,
        
        full_fit =
          full_ols,
        
        abnormality =
          abnormality,
        
        selected_idx =
          integer(0),
        
        affected_idx =
          affected_idx,
        
        runtime =
          full_mm_runtime
      )
    
    
    # ==========================================================================
    # 4.11 MIS
    #
    # Uses the FWL problem already prepared above.
    # ==========================================================================
    
    t_mis <-
      proc.time()[3L]
    
    
    mis <-
      run_mis_sensitivity(
        
        mod_full =
          mod,
        
        formula =
          formula,
        
        data =
          data,
        
        k =
          k,
        
        target_var =
          target_var,
        
        target_pos =
          target_pos,
        
        use_cpp =
          use_cpp,
        
        # --------------------------------------------------------------
        # Reuse full OLS
        # --------------------------------------------------------------
        
        full_fit =
          full_ols,
        
        # --------------------------------------------------------------
        # Reuse FWL projection
        # --------------------------------------------------------------
        
        problem =
          mis_problem
      )
    
    
    mis_runtime <-
      proc.time()[3L] -
      t_mis
    
    
    # ==========================================================================
    # 4.12 MIS -> OLS
    #
    # run_mis_sensitivity() has ALREADY fitted the selected deletion model.
    #
    # Therefore DO NOT run lm() again here.
    # ==========================================================================
    
    rows[[
      length(rows) + 1L
    ]] <-
      make_method_row05(
        
        base_meta =
          meta,
        
        model_state =
          state_name,
        
        diagnostic =
          "MIS",
        
        estimator =
          "OLS",
        
        fit =
          mis$deleted,
        
        full_fit =
          full_ols,
        
        abnormality =
          abnormality,
        
        selected_idx =
          mis$indices,
        
        affected_idx =
          affected_idx,
        
        runtime =
          mis_runtime
      )
    
    
    # ==========================================================================
    # 4.13 MIS -> MM
    #
    # This is retained because it answers:
    #
    #   Does MIS add anything beyond robust MM estimation?
    #
    # Only ONE MM deletion fit for MIS is done for each k.
    # ==========================================================================
    
    t_mis_mm <-
      proc.time()[3L]
    
    
    mis_mm <-
      fit_target_mm05(
        
        formula =
          formula,
        
        data =
          data,
        
        target_var =
          target_var,
        
        exclude_idx =
          mis$indices
      )
    
    
    mis_mm_runtime <-
      proc.time()[3L] -
      t_mis_mm
    
    
    rows[[
      length(rows) + 1L
    ]] <-
      make_method_row05(
        
        base_meta =
          meta,
        
        model_state =
          state_name,
        
        diagnostic =
          "MIS",
        
        estimator =
          "MM",
        
        fit =
          mis_mm,
        
        full_fit =
          full_ols,
        
        abnormality =
          abnormality,
        
        selected_idx =
          mis$indices,
        
        affected_idx =
          affected_idx,
        
        runtime =
          mis_mm_runtime
      )
    
    
    # ==========================================================================
    # 4.14 Classical top-k sets
    #
    # No diagnostics are recalculated here.
    #
    # We simply take:
    #
    #   first k ranked DFBETAS
    #   first k ranked Cook
    #   first k ranked leverage
    # ==========================================================================
    
    classical_sets <- list(
      
      DFBETAS =
        rank_dfb[
          seq_len(k)
        ],
      
      Cook =
        rank_cook[
          seq_len(k)
        ],
      
      Leverage =
        rank_lev[
          seq_len(k)
        ]
    )
    
    
    # ==========================================================================
    # 4.15 Classical diagnostics -> OLS
    #
    # Intentionally NOT fitting:
    #
    #   DFBETAS -> MM
    #   Cook     -> MM
    #   Leverage -> MM
    #
    # because these would add 12 lmrob() fits per model state across four k's.
    # ==========================================================================
    
    for (
      diagnostic_name
      in
      names(classical_sets)
    ) {
      
      idx <-
        classical_sets[[
          diagnostic_name
        ]]
      
      
      t_classical <-
        proc.time()[3L]
      
      
      fit_ols <-
        fit_target_ols05(
          
          formula =
            formula,
          
          data =
            data,
          
          target_var =
            target_var,
          
          exclude_idx =
            idx
        )
      
      
      classical_runtime <-
        proc.time()[3L] -
        t_classical
      
      
      rows[[
        length(rows) + 1L
      ]] <-
        make_method_row05(
          
          base_meta =
            meta,
          
          model_state =
            state_name,
          
          diagnostic =
            diagnostic_name,
          
          estimator =
            "OLS",
          
          fit =
            fit_ols,
          
          full_fit =
            full_ols,
          
          abnormality =
            abnormality,
          
          selected_idx =
            idx,
          
          affected_idx =
            affected_idx,
          
          runtime =
            classical_runtime
        )
    }
  }
  
  
  # ============================================================================
  # 4.16 Combine result rows
  # ============================================================================
  
  dplyr::bind_rows(
    rows
  )
}


# ==============================================================================
# 5. Run ONE Monte Carlo iteration
#
# A single iteration:
#
#   1. Generates ONE dataset
#   2. Runs intentionally misspecified model
#   3. Processes ALL k fractions
#   4. Optionally runs correctly respecified model
#
# k is NOT involved in generating the data.
# ==============================================================================

run_misspec_iter <- function(
    
  iter,
  
  n,
  
  x_type,
  
  error_type,
  
  scenario,
  
  severity_level,
  
  k_fractions = c(
    0.01,
    0.025,
    0.05,
    0.10
  ),
  
  seed_base =
    20260828L,
  
  use_cpp =
    TRUE,
  
  # --------------------------------------------------------------------------
  # DGP structural controls
  # --------------------------------------------------------------------------
  
  beta0 =
    0,
  
  beta1 =
    1,
  
  rho_xz =
    0.50,
  
  hetero_group_prop =
    0.25,
  
  break_fraction =
    0.75,
  
  threshold_quantile =
    0.75,
  
  # --------------------------------------------------------------------------
  # Distribution controls
  # --------------------------------------------------------------------------
  
  mix_prop =
    0.10,
  
  skew_t_df =
    5,
  
  gpd_shape =
    0.25,
  
  pareto_shape =
    3
  
) {
  
  # ============================================================================
  # 5.1 Stable Monte Carlo seed
  # ============================================================================
  
  set.seed(
    as.integer(
      seed_base +
        iter
    )
  )
  
  
  # ============================================================================
  # 5.2 Generate ONE dataset
  #
  # Notice:
  #   k_fraction does NOT appear here.
  # ============================================================================
  
  dat <-
    generate_misspec_data(
      
      n =
        n,
      
      x_type =
        x_type,
      
      error_type =
        error_type,
      
      scenario =
        scenario,
      
      severity_level =
        severity_level,
      
      beta0 =
        beta0,
      
      beta1 =
        beta1,
      
      rho_xz =
        rho_xz,
      
      hetero_group_prop =
        hetero_group_prop,
      
      break_fraction =
        break_fraction,
      
      threshold_quantile =
        threshold_quantile,
      
      mix_prop =
        mix_prop,
      
      skew_t_df =
        skew_t_df,
      
      gpd_shape =
        gpd_shape,
      
      pareto_shape =
        pareto_shape
    )
  
  
  # ============================================================================
  # 5.3 Metadata that does NOT depend on k
  # ============================================================================
  
  base_meta <- list(
    
    iter =
      iter,
    
    n =
      n,
    
    true_beta =
      dat$true_beta,
    
    draw_seed =
      as.integer(
        seed_base + iter
      ),
    
    x_type =
      x_type,
    
    error_type =
      error_type,
    
    scenario =
      scenario,
    
    severity_level =
      severity_level,
    
    severity_target =
      dat$severity_target,
    
    severity_realized =
      dat$severity_realized,
    
    raw_misspec_coef =
      dat$raw_misspec_coef
  )
  
  
  # ============================================================================
  # 5.4 Intentionally misspecified model
  #
  # ALL k fractions operate on the SAME data and SAME full model.
  # ============================================================================
  
  wrong <-
    run_state05(
      
      data =
        dat$data,
      
      formula =
        dat$wrong_formula,
      
      state_name =
        "wrong",
      
      k_fractions =
        k_fractions,
      
      affected_idx =
        dat$affected_idx,
      
      base_meta =
        base_meta,
      
      target_var =
        "x",
      
      use_cpp =
        use_cpp
    )
  
  
  # ============================================================================
  # 5.5 Correct baseline scenario
  #
  # Under scenario == "correct":
  #
  #   wrong_formula == correct_formula
  #
  # so running the model again would be pure duplication.
  # ============================================================================
  
  if (
    scenario ==
    "correct"
  ) {
    
    return(
      wrong
    )
  }
  
  
  # ============================================================================
  # 5.6 Endogeneity and heteroskedasticity
  #
  # These are not repaired simply by adding an omitted term to the regression
  # formula, so do not create a misleading "correct specification" state.
  # ============================================================================
  
  if (
    scenario %in%
    c(
      "endogeneity",
      "heteroskedastic"
    )
  ) {
    
    return(
      wrong
    )
  }
  
  
  # ============================================================================
  # 5.7 Correct structural specification
  #
  # Applies to:
  #
  #   OVB
  #   nonlinearity
  #   heterogeneous slope
  #   structural break
  #   threshold
  #   missing interaction
  # ============================================================================
  
  correct <-
    run_state05(
      
      data =
        dat$data,
      
      formula =
        dat$correct_formula,
      
      state_name =
        "correct",
      
      k_fractions =
        k_fractions,
      
      affected_idx =
        dat$affected_idx,
      
      base_meta =
        base_meta,
      
      target_var =
        "x",
      
      use_cpp =
        use_cpp
    )
  
  
  # ============================================================================
  # 5.8 Return both model states
  # ============================================================================
  
  dplyr::bind_rows(
    wrong,
    correct
  )
}