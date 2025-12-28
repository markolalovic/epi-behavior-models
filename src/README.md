# Source directory
This directory contains implementation of the models, the ABC-SMC inference pipeline, and post-processing scripts to create data for plotting figures.

## Core:
- `config.py`: configuration for locations, priors, and hyper-parameters
- `models.py`: model definitions and model runner solving system of ODEs
- `inference.py`: distance function definition and ABC-SMC runner
- `utils.py`: some math and data helpers
- `plotting.py`: plotting helpers for KDEs and fits figures

## Runners:
- `run_inference.py`: estimates parameter posteriors for all 30 locations
- `run_selection.py`: performs pairwise model selection (Bayes factors)
- `run_synthetic_experiments.py`: structural bias experiments on synthetic data
- `run_example.py`: quick parameter recovery examples for main README.md

## Pre-processing scripts:
- `process_mortality_data.py`: to process the raw JHU data.

## Post-processing scripts:
To process results from `results/` and create CSVs in `data/plotting/`:
- `create_fits_data.py`: medians and stats for fits figures and uncertainty bands
- `create_r0_data.py, create_final_epi_size_data.py`: stats for boxplot figures
- `median_nssr_distances.py`: weighted medians of distances as metrics for comparison tables
- `create_synthetic_*_data.py`: stats and fits for synthetic experiments
- `create_beta_eff_data.py, create_re_data.py`: extra time dependent transmission rates and reproduction numbers for supplement figures
