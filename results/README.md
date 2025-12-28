# Results directory
This directory stores the output of the ABC-SMC inference.

## Examples
The estimated posteriors for 2 examples are stored in:
- `results/example_baseline/baseline_posterior_1.csv`
- `results/example_behavior/behavior_mixed_posterior_1.csv`

From running `src/run_example.py`:
```bash
python3 run_example.py --model baseline
python3 run_example.py --model behavior
```

## Summary results
- Weighted median NSSR distances: `results/model_comparison/median_nssr_distances.csv`
  * See [`src/median_nssr_distances.py`](../src/median_nssr_distances.py) for details.

- Model selection results summaries: `results/model_selection/baseline_vs_{variant}_summary_1.csv`
  * See [`src/run_selection.py`](../src/run_selection.py) for details.

- Tables for main and supp: `results/table_comparison_baseline_{variant}.tex`
  * See [`R/generate_tables.R`](../R/generate_tables.R) for details.

## Parameter inference
- Estimated posteriors as samples and weights: `results/{location}/{model}_posterior_{RESULT}.csv`

- Estimated posteriors over generations (full pyABC history object): `results/{location}/{model}_history_{RESULT}.pkl`

Where:
- `location` is US state abbreviation from `LOCATIONS` from `src/config.py`.

- `model` refers to baseline or one of behavioral variants: `behavior_exp`, `behavior_mixed`, `behavior_rational` defined in `src/models.py`

- `RESULT` is the encoding from `src/config.py` (e.g., 1 is the first run, 2 the second run etc.)

These outputs are not tracked by version control to keep repository lightweight; see `.gitignore`.

To recreate them run [src/run_inference.py](../src/run_inference.py).

This runs the loop over 30 locations and 4 models and stores the outputs to `results/` directory.

## Model selection
- As for parameter inference, posterior samples, weights and history objects are saved to: `results/model_selection/`

- Convergence plots of posterior model probabilities P(M | D, t) over generations: `results/model_selection/{location}_{behavior_name}_trajectory_{RESULT}.png`

These outputs are not tracked by version control to keep repository lightweight; see `.gitignore`.

To recreate them run [src/run_selection.py](../src/run_selection.py).







