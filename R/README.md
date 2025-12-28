# R visualization directory
This directory contains the R scripts to generate the figures and tables in the main text and Supplement.

## Main:
- `generate_tables.R`: merges results into: 
  * `results/table_comparison_baseline_exp.tex`
  * `results/table_comparison_baseline_mixed.tex`
  * `results/table_comparison_baseline_rational.tex`

- `plot_fits_data.R`: model fits grid figures comparing Baseline and behavioral variants for all 30 locations

- `plot_r0_and_final_size_boxplots.R`: 2-panel boxplots comparing $R_0$ and final epidemic sizes

- `plot_synthetic_results.R`: 3-panel recovery figures for structural bias experiments

## Supplement:
- `plot_bad_locations.R`: mortality data for excluded states

- `plot_beta_eff_all.R`: time dependent effective transmission rates $\beta(t)$

- `plot_re_eff_all.R`: time dependent effective reproduction numbers $R_e(t)$

