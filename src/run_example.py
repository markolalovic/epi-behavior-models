#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/run_example.py
For example model:
- Generates synthetic data daily deaths using config.EXAMPLE_PARAMS
- Runs ABC-SMC to estimate posteriors
- Plots posteriors and model fit to synthetic data

For Baseline model (SEIRD) run:
    python3 run_example.py --model baseline

For Behavior model (Mixed) run:
    python3 run_example.py --model behavior

Outputs:
- figures/examples/fit_{label}.png"

Where label is either `baseline` or `behavior`.

Additional results and diagnostics are stored to:
- results/example_{label}/{model_name}_posterior_1.csv
- data/local-db/example_{label}_{model_name}_{RESULT}.db
- results/example_{label}/{model_name}_history_1.pkl
"""

import os
import argparse
import numpy as np
import pandas as pd

from config import (
    MODELS, FIXED_PARAMS, PRIORS, START_DATE, END_DATE, 
    EXAMPLE_N, EXAMPLE_PARAMS, RESULT
)
from run_synthetic_experiments import generate_synthetic_data
from inference import run_abc
from plotting import setup_plotting, plot_posteriors, plot_fits

if __name__ == "__main__":
    # parsing model argument
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["baseline", "behavior"], required=True)
    args = parser.parse_args()

    # configure model dictionary
    if args.model == "behavior":
        model_dict = next(m for m in MODELS if m["name"] == "behavior_mixed")
        truth = EXAMPLE_PARAMS
        beta_form = "mixed"
    else:
        model_dict = next(m for m in MODELS if m["name"] == "baseline")
        truth = {k: v for k, v in EXAMPLE_PARAMS.items() if k != "zeta"}
        beta_form = "constant"

    print(f"\n{'='*60}\nRunning for {args.model}\n{'='*60}")
    
    print(f"Generating synthetic data...")
    obs_length = (pd.to_datetime(END_DATE) - pd.to_datetime(START_DATE)).days + 1
    fixed = FIXED_PARAMS.copy()
    fixed["N"] = EXAMPLE_N
    
    np.random.seed(12345)
    obs_deaths = generate_synthetic_data(
        beta_form=beta_form,
        base_params=truth,
        zeta_value=truth.get("zeta", 0.0),
        fixed_params=fixed,
        obs_length=obs_length
    )

    result_location = f"example_{args.model}"
    post_path = f"../results/{result_location}/{model_dict['name']}_posterior_{RESULT}.csv"
    if not os.path.exists(post_path):
        print(f"Estimating model parameters...")
        history = run_abc(
            fixed_params=fixed,
            obs_deaths=obs_deaths,
            base_prior=PRIORS,
            location=result_location,
            model_dict=model_dict
        )
        df_post, w_post = history.get_distribution(m=0, t=history.max_t)
    else:
        # posterior results are already there 
        print(f"Loading existing results from: {post_path}")
        df_full = pd.read_csv(post_path)
        w_post = df_full["weight"].values
        df_post = df_full.drop(columns=["weight"])

    # plot results
    setup_plotting(font_size=15)
    plot_posteriors(df_post, w_post, truth, args.model)
    plot_fits(df_post, w_post, model_dict, obs_deaths, args.model)
    print("\n Done.")
