# MIS vs. Classical Diagnostics — Phase 1

Benchmarks classical leave-one-out (LOO) regression diagnostics (Cook's Distance,
DFBETAS, Leverage) against exact Most Influential Sets (MIS) detection under two
failure modes: **masking** and **sign-flipping**.

---

## Repository Structure

| Path | Role |
|---|---|
| `R/diagnostics_classical.R` | LOO baseline toolkit |
| `script/01_classical_vs_mis_detection.R` | Simulation script |
| `output/fig1_classical_vs_mis_detection.png` | Output figure |

---

## R Module: `diagnostics_classical.R`

Standardised interface for computing LOO diagnostics, designed for direct
comparison against MIS output.

| Function | Description |
|---|---|
| `get_leverage(model)` | Hat matrix diagonal |
| `get_cooks_d(model)` | Cook's Distance for all observations |
| `get_dfbetas(model, target_var)` | DFBETAS for a specific covariate |
| `get_all_classical(model, target_var)` | All metrics in a single `data.frame` |
| `get_classical_set(model, target_var, k, metric)` | Top-k set or threshold-based set ($4/n$, $2p/n$, $2/\sqrt{n}$) |

### Usage

```r
source("R/diagnostics_classical.R")
fit <- lm(mpg ~ wt + hp, data = mtcars)

# All diagnostics in one table
metrics_df <- get_all_classical(fit, target_var = "wt")

# Top 3 by Cook's Distance
top3 <- get_classical_set(fit, target_var = "wt", k = 3, metric = "cooks_d")

# Threshold-based DFBETAS flagging
flagged <- get_classical_set(fit, target_var = "wt", k = NULL, metric = "dfbetas_target")
```

---

## Script: `01_classical_vs_mis_detection.R`

Constructs a synthetic DGP with a true positive slope ($Y = 1.5X + \varepsilon$),
then injects two pathological structures:

- **Masked pairs** — tightly clustered outliers whose marginal LOO influence is
  near zero, blinding Cook's D, DFBETAS, and Leverage.
- **Sign-flippers** — a high-leverage cluster that inverts the estimated slope
  from positive to negative when classical tools fail to flag it.

Runs all four detection methods (top $k = 10$) and writes a 2×2 figure comparing
flagged points, MIS-only flags, tool agreement, and post-removal fitted lines.