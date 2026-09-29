# Most Influential Sets: Sensitivity in Statistical Models and Applications

This repository contains the code, simulation results, empirical outputs, and LaTeX source for my Master of Applied Econometrics thesis at Monash University.

The thesis studies finite-sample sensitivity of regression coefficients using fixed-\(k\) Most Influential Sets (MIS).

## Replication

The main simulation workflows are in `scripts/`, reusable functions are in `R/`, and generated results are stored in `output/`.

The code includes:

- exact fixed-\(k\) MIS search for a prespecified regression coefficient;
- optional Rcpp acceleration;
- comparisons with leverage, Cook's distance, DFBETAS, MM, and LTS;
- simulations for contamination and model misspecification;
- full-refit validation after MIS selection;
- parallel computation and resumable checkpoint files;
- scripts for generating thesis tables and figures.

Some simulations are computationally intensive and may take many hours.

For reference, development runs using 12 parallel workers included:

| Workflow | Approximate wall time |
|---|---:|
| `05a` null calibration | 2.5 hours |
| `05a` evaluation | 20 hours |
| `05c` correct-model calibration | 2.5 hours |
| `05c` with classical diagnostics | 6 hours |

Actual runtime depends on CPU allocation, memory, worker count, and whether the Rcpp implementation is enabled.

The simulation scripts should generally be run from the `scripts/` directory.

## Repository Structure

```text
R/             Reusable statistical and simulation functions
src/           Optional Rcpp acceleration
scripts/       Simulation and output-generation scripts
output/        Generated simulation results, tables, and figures
application/   Empirical application outputs
paper/         LaTeX thesis source