# Most Influential Sets: Sensitivity in Statistical Models and Applications

This repository contains the code, simulation results, empirical outputs, and LaTeX source for my ETC5860 Honours Research Project, completed as part of the Master of Applied Econometrics at Monash University.

The project studies finite-sample sensitivity of regression coefficients using fixed-\(k\) Most Influential Sets (MIS).

## Computational Contribution

The repository develops a reproducible computational framework for applying and evaluating fixed-\(k\) MIS. It includes:

- an R implementation of coefficient-targeted fixed-\(k\) MIS with Frisch-Waugh-Lovell residualisation and optional Rcpp acceleration;
- simulation frameworks comparing MIS with classical influence diagnostics and robust estimators under contamination and model misspecification;
- full-refit validation, parallel computation, checkpointing, and automated generation of thesis tables and figures.

The underlying fixed-\(k\) MIS formulation and Dinkelbach optimisation are based on existing methodological work.

## Replication

The main simulation workflows are in `scripts/`, reusable functions are in `R/`, and generated results are stored in `output/`.

Some simulations are computationally intensive and may take many hours. For reference, development runs used 12 parallel CPU workers on an Intel Xeon Gold 6548Y+ system (up to 4.1 GHz) with 1 TiB of RAM; the `05a` evaluation took approximately 20 hours.

The simulation scripts should generally be run from the `scripts/` directory.

## Repository Structure

```text
R/             Reusable statistical and simulation functions
src/           Optional Rcpp acceleration
scripts/       Simulation and output-generation scripts
output/        Generated simulation results, tables, and figures
application/   Empirical application outputs
paper/         LaTeX thesis source