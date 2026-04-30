#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/run_synthetic_experiments.py
Runs synthetic experiments with different ground truths:
  - Behavior (Mixed)
  - Behavior (Exponential)
  - Behavior (Rational)

For each beta_form in {"mixed", "exp", "rational"}:
  For each zeta value in a grid of zeta values:
  - Generates synthetic smoothed mortality series

Where:
  - Other parameter values are weighted medians from MA:
    - behavior_{beta_form} posterior
  - Gaussian noise is added to simulated extracted daily deaths series
  - Data is aligned to chosen dates and trailing MA(7) is applied,
  - same as the for the real data preprocessing.

For each synthetic dataset, estimate parameter posteriors for: 
  - Baseline and all behavior variants.

Outputs:
- Synthetic datasets with columns: `date, y_zeta_0.001, ..., y_zeta_0.02`:
  * `data/synthetic/synthetic_mortality_ground_truth_{beta_form}.csv`

- Histories and estimated posteriors are saved to:
  * `results/synthetic_ground_truth_{beta_form}/synthetic_zeta{k}/`

- As files:
  * `baseline_history_{RESULT}.pkl`
  * `baseline_posterior_{RESULT}.csv`
  * `behavior_{beta_form}_history_{RESULT}.pkl`
  * `behavior_{beta_form}_posterior_{RESULT}.csv`
"""

import os
import numpy as np
import pandas as pd

from config import (
    MODELS, FIXED_PARAMS, POPULATION_SIZE, 
    PRIORS, SYNTHETIC_CONFIG, START_DATE, END_DATE
)

from models import run_seird_model
from inference import run_abc
from utils import trailing_ma7, variant_name, ground_truth_values

# synthetic experiment configuration
gt_location = SYNTHETIC_CONFIG["ground_truth_location"]
zeta_grid = SYNTHETIC_CONFIG["zeta_grid"]
noise_SD = SYNTHETIC_CONFIG["noise_sd"]

def generate_synthetic_data(beta_form, base_params, zeta_value, 
                            fixed_params, obs_length, smoothing_window_k=7):
    """
    Generates synthetic data with added noise and smoothing. 
    Set `beta_form = "constant"` to generate data under the baseline model.
    """
    if beta_form == "constant":
        model_name = "baseline"
        display_name = "Baseline"
        is_behavioral = False
    else:
        model_name = variant_name(beta_form)
        display_name = f"Behavioral ({beta_form})"
        is_behavioral = True
    
    model_dict = {
        "name": model_name,
        "display_name": display_name,
        "is_behavioral": is_behavioral,
        "beta_form": beta_form,
    }

    params = dict(base_params)
    if not is_behavioral:
        params["zeta"] = 0.0
    else:
        params["zeta"] = zeta_value
    
    sim = run_seird_model(
        params=params,
        fixed_params=fixed_params,
        obs_length=obs_length,
        model_dict=model_dict
    )

    # extract raw daily deaths
    traj = sim["full_trajectory"]
    daily_raw = np.diff(traj[:, 4]) # length = obs_length + pad = 123 + 6

    # add Gaussian noise
    sigma = noise_SD * max(1.0, np.max(daily_raw))
    noise = np.random.normal(0.0, sigma, size=len(daily_raw))
    noisy_raw = daily_raw + noise

    # remove any negative values
    noisy_raw = np.maximum(noisy_raw, 0.0)

    # apply trailing MA(7) smoothing
    noisy_smoothed = trailing_ma7(noisy_raw)
    
    # align to target window by dropping the leading zeros
    pad = smoothing_window_k - 1
    final_obs = noisy_smoothed[pad:pad + obs_length]

    return final_obs

def generate_synthetic_data_all_zeta_values(beta_form):
    dates = pd.date_range(
        pd.to_datetime(START_DATE),
        pd.to_datetime(END_DATE),
        freq="D")
    obs_length = len(dates)

    fixed_params = FIXED_PARAMS.copy()
    fixed_params["N"] = POPULATION_SIZE[gt_location]

    base_params = ground_truth_values(beta_form)

    data = {"date": dates}
    for z in zeta_grid:
        y = generate_synthetic_data(beta_form, base_params, z, fixed_params, obs_length=obs_length)
        col = f"y_zeta_{z:g}"
        data[col] = y

    wide_df = pd.DataFrame(data)
    return dates, wide_df, {"base_params": base_params, "fixed_params": fixed_params}

def fit_to_synthetic_data(synthetic_data, location_path, models_to_fit):
    """
    Fits models_to_fit to synthetic_data:
      - synthetic_data: smoothed simulated daily deaths
      - models_to_fit from config.MODELS
      - location_path: `results/synthetic_ground_truth_{beta_form}`
    """
    fixed_params = FIXED_PARAMS.copy()
    fixed_params["N"] = POPULATION_SIZE[gt_location]
    obs_deaths = np.asarray(synthetic_data, dtype=float)

    for model_dict in models_to_fit:
        print(f"Fitting {model_dict['name']} to synthetic_data...")
        prior = PRIORS.copy()
        run_abc(
            fixed_params=fixed_params,
            obs_deaths=obs_deaths,
            base_prior=prior,
            location=location_path,
            model_dict=model_dict,
        )

if __name__ == "__main__":
    print("Running synthetic experiments for all behavior ground truths.")
    np.random.seed(12345)

    for beta_form in ["exp", "mixed", "rational"]:
        print(f"\nRunning for ground truth: {beta_form}...")

        # Generate synthetic dataset
        print("\nGenerating synthetic dataset...")
        dates, wide_df, meta = generate_synthetic_data_all_zeta_values(beta_form)
        data_out_wide = f"../data/synthetic_mortality_ground_truth_{beta_form}.csv"
        wide_df.to_csv(data_out_wide, index=False)
        print(f"Saved synthetic synthetic dataset to: {data_out_wide}")

        print(f"Ground truth weighted-medians for {gt_location} with behavior variant {beta_form}:")
        for k, v in meta["base_params"].items():
            print(f"    {k}: {v:.6g}")
        print(f"N: {meta['fixed_params']['N']}")

        # fit all models to this synthetic dataset
        print("\nFitting all models to synthetic dataset...")
        for i, z in enumerate(zeta_grid, start=1):
            col = f"y_zeta_{z:g}"
            synthetic_data = wide_df[col]

            print(f"Fitting to synthetic_data generated with zeta = {z}...")
            zeta_tag = f"synthetic_zeta{i}"
            # we need to pass to run_abc the path where results should be saved
            # results/synthetic_ground_truth_{beta_form}/sythetic_zeta{number}
            location_path = os.path.join(f"synthetic_ground_truth_{beta_form}", zeta_tag)
            fit_to_synthetic_data(synthetic_data, location_path, MODELS)

        print(f"\nDone for ground truth: {beta_form}.")

    print("\nDone for all ground truths.")

