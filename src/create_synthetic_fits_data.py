#!/usr/bin/env python3
# src/create_synthetic_fits_data.py
#
# Generates the quantiles for synthetic fits figures Panels A.
# Calculates: 
# - medians and 90% credible intervals from 1000 simulations per synthetic experiment.
#
# For each baseline, behavior variant(beta_form) pair,
# where beta_form in ["mixed", "exp", "rational"]:
# - Load synthetic smoothed synthetic observations for zeta values:
#     `data/synthetic_mortality_ground_truth_{beta_form}.csv`
#
# - For each zeta value:
#    - load weighted posterior
#    - sample 1000 parameter sets (with weights)
#    - simulate trajectories
#    - compute median and 90% credible interval
#    - save columns:
#     date, synthetic_id, zeta_value, beta_form, model_name, 
#     lower_90, median, upper_90, synth observation
#
# Outputs:
# - data/plotting/fits_data_synthetic_ground_truth_mixed.csv
# - data/plotting/fits_data_synthetic_ground_truth_exp.csv
# - data/plotting/fits_data_synthetic_ground_truth_rational.csv
#

import os
import numpy as np
import pandas as pd
from config import MODELS, POPULATION_SIZE, FIXED_PARAMS, RESULT, SYNTHETIC_CONFIG
from models import run_seird_model

def generate_predictive_trajectories(posterior_df, model_dict, obs_length, 
                                     n_samples=1000):
    """Generates ensemble of simulated trajectories by sampling 
    from estimated posterior using weights."""
    sampled = posterior_df.sample(
        n=n_samples, 
        weights="weight", 
        replace=True).to_dict("records")
    
    gt_location = SYNTHETIC_CONFIG["ground_truth_location"]
    fixed_params = FIXED_PARAMS.copy()
    fixed_params["N"] = POPULATION_SIZE[gt_location]

    sims = []
    for params in sampled:
        out = run_seird_model(params, fixed_params, obs_length, model_dict)
        sims.append(out["data"])
    
    return np.array(sims, dtype=float)

def generate_fits_data_for_ground_truth(beta_form):
    """Generates fits data for specific ground truth variant mixed, exp, or rational."""
    print(f"\nCreating synthetic fits data for ground truth: {beta_form}")

    # data paths
    syn_data_path = f"../data/synthetic/synthetic_mortality_ground_truth_{beta_form}.csv"
    out_path = f"../data/plotting/fits_data_synthetic_ground_truth_{beta_form}.csv"
    base_res_dir = f"../results/synthetic_ground_truth_{beta_form}"

    # load synthetic observations
    df_wide = pd.read_csv(syn_data_path, parse_dates=["date"]).sort_values("date")
    dates = df_wide["date"].values
    
    # identify baseline and specific behavior variant used for truth
    baseline_model = next(m for m in MODELS if m["name"] == "baseline")
    behavior_model = next(m for m in MODELS if m.get("beta_form") == beta_form)
    models_to_process = [baseline_model, behavior_model]

    # generate blocks of data for all 6 zeta values in the synthetic grid
    all_blocks = []

    # set columns as `y_zeta_0.001, y_zeta_0.005, ...`
    series_cols = [c for c in df_wide.columns if c.startswith("y_zeta_")]
    series_cols_sorted = sorted(series_cols, key=lambda x: float(x.split('_')[-1]))

    for i, col in enumerate(series_cols_sorted, start=1):
        y_obs = df_wide[col].to_numpy()
        zeta_val = float(col.split('_')[-1])
        zeta_tag = f"synthetic_zeta{i}"

        for m_dict in models_to_process:
            m_name = m_dict["name"]
            post_path = f"{base_res_dir}/{zeta_tag}/{m_name}_posterior_{RESULT}.csv"

            if not os.path.exists(post_path):
                continue

            post_df = pd.read_csv(post_path)
            
            # generate 1000 trajectories
            trajs = generate_predictive_trajectories(post_df, m_dict, len(y_obs))
            
            # get the quantiles for the bands and median
            q_low, q_med, q_hi = np.quantile(trajs, [0.05, 0.5, 0.95], axis=0)

            # construct the data block
            block = pd.DataFrame({
                "date": dates,
                "zeta_true": zeta_val,
                "model_name": m_dict.get("display_name", m_name),
                "lower_90": q_low,
                "median": q_med,
                "upper_90": q_hi,
                "obs": y_obs
            })
            all_blocks.append(block)
        
        print(f"Processed {zeta_tag} zeta = {zeta_val}.")

    # concatenate all blocks and save it
    final_df = pd.concat(all_blocks, ignore_index=True)
    final_df.to_csv(out_path, index=False)
    print(f"Saved synthetic fits data to {out_path}")

if __name__ == "__main__":
    for variant in ["mixed", "exp", "rational"]:
        generate_fits_data_for_ground_truth(variant)
