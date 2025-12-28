#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/create_re_data.py
Calculates the data for plotting the effective reproduction number Re(t)

Outputs:
- data/plotting/re_effective_data.csv
"""

import numpy as np
import pandas as pd

from config import MODELS, LOCATIONS, POPULATION_SIZE, FIXED_PARAMS, RESULT
from models import run_seird_model

def re_series(model_dict, R0_base, zeta, delta, S_series, I_series, N):
    """Calculates Re(t) = R0(t) * S(t)/N given behavioral form."""
    if not model_dict["is_behavioral"]:
        R0_t = np.full_like(I_series, R0_base)
    else:
        # M(t) = delta * I(t)
        M = delta * I_series
        arg = zeta * M
        
        # clip for stability as in models.py
        arg_exp = np.clip(arg, -50.0, 50.0)
        arg_den = np.clip(arg, -0.99, 50.0)

        beta_form = model_dict["beta_form"]
        if beta_form == "exp":
            R0_t = R0_base * np.exp(-arg_exp)
        elif beta_form == "rational":
            R0_t = R0_base / (1.0 + arg_den)
        elif beta_form == "mixed":
            R0_t = R0_base * np.exp(-arg_exp) / (1.0 + arg_den)
        else:
            R0_t = np.full_like(I_series, R0_base)

    return R0_t * (S_series / N)

def generate_re_trajectories(post_df, model_dict, loc, obs_length, n_samples=1000):
    """Generates ensemble of Re(t) trajectories by sampling from estimated posterior."""
    samples = post_df.sample(n=n_samples, weights="weight", replace=True).to_dict("records")
    
    fixed = FIXED_PARAMS.copy()
    fixed["N"] = POPULATION_SIZE[loc]
    N = fixed["N"]
    
    re_curves = []
    for row in samples:
        # run simulation to get S(t) and I(t)
        sim = run_seird_model(row, fixed, obs_length, model_dict)
        traj = sim["full_trajectory"] # length obs_length + pad + 1
        
        # extract compartments skipping leading pad days
        pad = 6
        S_aligned = traj[pad : pad + obs_length, 0] # S(t) at index 0
        I_aligned = traj[pad : pad + obs_length, 2] # I(t) at index 2
        
        # params
        delta_nat = np.exp(row["delta"])
        R0_base = row["R0"]
        zeta = row.get("zeta", 0.0)
        
        # calculate Re(t)
        curve = re_series(model_dict, R0_base, zeta, delta_nat, S_aligned, I_aligned, N)
        re_curves.append(curve)
        
    return np.array(re_curves)

def create_re_data():
    print("Generating effective reproduction number Re plotting data...")
    
    # data paths
    df_obs = pd.read_csv("../data/processed/smoothed_mortality.csv", parse_dates=["date"])
    dates, obs_length = df_obs["date"].values, len(df_obs)

    all_results = []
    for loc in LOCATIONS:
        print(f"Processing {loc}...")
        for model in MODELS:
            post_path = f"../results/{loc}/{model['name']}_posterior_{RESULT}.csv"
            post_df = pd.read_csv(post_path)
            
            # get ensemble and compute 90% uncertainty intervals
            re_matrix = generate_re_trajectories(post_df, model, loc, obs_length)
            q_low, q_med, q_hi = np.quantile(re_matrix, [0.05, 0.5, 0.95], axis=0)
            
            all_results.append(pd.DataFrame({
                "date": dates,
                "location": loc,
                "model_name": model["display_name"],
                "median": q_med,
                "lower": q_low,
                "upper": q_hi
            }))

    final_df = pd.concat(all_results, ignore_index=True)
    out_path = "../data/plotting/re_effective_data.csv"
    
    final_df.to_csv(out_path, index=False)
    print(f"Saved to {out_path}")

if __name__ == "__main__":
    create_re_data()
