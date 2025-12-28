#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/create_beta_eff_data.py
Calculates the data for plotting the effective transmission rate beta(t).

Outputs:
- data/plotting/beta_eff_data.csv

"""

import os
import numpy as np
import pandas as pd
from config import MODELS, LOCATIONS, POPULATION_SIZE, FIXED_PARAMS, RESULT
from models import run_seird_model
from utils import weighted_quantile

def beta_eff_series(beta_form, beta0, zeta, delta, I_series):
    """Calculates beta(t) series given beta_form."""
    if beta_form == "constant":
        return np.full_like(I_series, beta0)
    
    # M(t) = delta * I(t)
    M = delta * I_series
    arg = zeta * M
    
    # for stability clip as in models.py
    arg_exp = np.clip(arg, -50.0, 50.0)
    arg_den = np.clip(arg, -0.99, 50.0)

    if beta_form == "exp":
        return beta0 * np.exp(-arg_exp)
    elif beta_form == "rational":
        return beta0 / (1.0 + arg_den)
    elif beta_form == "mixed":
        return beta0 * np.exp(-arg_exp) / (1.0 + arg_den)
    return np.full_like(I_series, beta0)

def generate_beta_trajectories(post_df, model_dict, loc, obs_length, n_samples=1000):
    """Simulates bundle of trajectories and calculates beta(t) for each."""
    samples = post_df.sample(n=n_samples, weights="weight", replace=True).to_dict("records")
    
    fixed = FIXED_PARAMS.copy()
    fixed["N"] = POPULATION_SIZE[loc]
    gamma = fixed["gamma"]
    
    beta_curves = []
    for row in samples:
        # run full SEIRD model to get I(t)
        sim = run_seird_model(row, fixed, obs_length, model_dict)
        traj = sim["full_trajectory"] # with padding of length obs_length + 6 + 1
        
        # align I(t) to target window
        # run_seird_model simulates obs_length + pad days. 
        # drop the first 'pad' indices to match the smoothed data
        pad = 6 
        I_aligned = traj[pad : pad + obs_length, 2] # I(t) is at index 2
        
        # reconstruct beta parameters after reparametrization
        delta_nat = np.exp(row["delta"])
        beta0 = row["R0"] * (gamma + delta_nat)
        zeta = row.get("zeta", 0.0)
        
        # calculate beta(t)
        curve = beta_eff_series(model_dict["beta_form"], beta0, zeta, delta_nat, I_aligned)
        beta_curves.append(curve)
    
    return np.array(beta_curves)

def create_beta_eff_data():
    print("Generating effective transmission rate beta(t) data...")
    
    # data paths
    data_path = "../data/processed/smoothed_mortality.csv"
    df_obs = pd.read_csv(data_path, parse_dates=["date"]).sort_values("date")
    dates, obs_length = df_obs["date"].values, len(df_obs)

    all_results = []
    for loc in LOCATIONS:
        print(f"Processing {loc}...")
        for model in MODELS:
            post_path = f"../results/{loc}/{model['name']}_posterior_{RESULT}.csv"
            if not os.path.exists(post_path): continue
            
            post_df = pd.read_csv(post_path)
            
            # get ensemble of curves
            beta_matrix = generate_beta_trajectories(post_df, model, loc, obs_length)
            
            # get pointwise 5th, 50th, 95th quantiles
            q_low, q_med, q_hi = np.quantile(beta_matrix, [0.05, 0.5, 0.95], axis=0)
            
            all_results.append(pd.DataFrame({
                "date": dates,
                "location": loc,
                "model_name": model["display_name"],
                "median": q_med,
                "lower": q_low,
                "upper": q_hi
            }))

    final_df = pd.concat(all_results, ignore_index=True)
    out_path = "../data/plotting/beta_eff_data.csv"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    final_df.to_csv(out_path, index=False)
    print(f"Saved to {out_path}")

if __name__ == "__main__":
    create_beta_eff_data()
