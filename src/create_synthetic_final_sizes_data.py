#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/create_synthetic_final_sizes_data.py
Calculated data for synthetic experiments Panel C figures:
- Calculates extinction-based final epidemic size:
    final_size = 1 - S(t_end) / N
  where t_end is the first day such that I(t_end) < 1.

For each ground truth beta_form in ["mixed", "exp", "rational"]:
- Loads synthetic mortality:
  `data/synthetic/synthetic_mortality_ground_truth_{beta_form}.csv
- And results synthetic_zeta{idx} results for all models:
  `results/synthetic_ground_truth_{beta_form}/synthetic_zeta{idx}/{model}_posterior_{RESULT}.csv
- Draws final size samples and computes boxplot stats.

Outputs:
- data/plotting/final_size_boxplot_data_synthetic_ground_truth_mixed.csv
- data/plotting/final_size_boxplot_data_synthetic_ground_truth_exp.csv
- data/plotting/final_size_boxplot_data_synthetic_ground_truth_rational.csv

"""

import os
import numpy as np
import pandas as pd
from config import MODELS, RESULT, POPULATION_SIZE, FIXED_PARAMS, SYNTHETIC_CONFIG
from models import compute_final_size_until_extinction

def calculate_extinction_stats(posterior_df, model_dict, n_samples=1000):
    """Calculates extinction-based final sizes by sampling from estimated posteriors."""
    # draw 1000 parameter sets using weights
    draws = posterior_df.sample(n=n_samples, weights="weight", replace=True).to_dict("records")
    
    fixed_params = FIXED_PARAMS.copy()
    gt_location = SYNTHETIC_CONFIG["ground_truth_location"]
    fixed_params["N"] = POPULATION_SIZE[gt_location]
    
    final_sizes = []
    for params in draws:
        _, _, fs = compute_final_size_until_extinction(
            params=params,
            fixed_params=fixed_params,
            model_dict=model_dict,
            total_days=1825 # ~ 5 years  
        )
        final_sizes.append(fs)
    
    sample = np.sort(np.asarray(final_sizes, dtype=float))
    
    # calculate boxplot stats
    q1, med, q3 = np.quantile(sample, [0.25, 0.5, 0.75])
    iqr = q3 - q1
    low_bound, high_bound = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    
    # whiskers are min / max values !within the prior bounds!
    xin = sample[(sample >= low_bound) & (sample <= high_bound)]
    return {
        "ymin": float(np.min(xin)),
        "lower": float(q1),
        "middle": float(med),
        "upper": float(q3),
        "ymax": float(np.max(xin))
    }

def generate_final_size_data_for_ground_truth(beta_form):
    """Processes final size samples for given behavioral ground truth."""
    print(f"\nCreating synthetic final size data for ground truth: {beta_form}")

    # data paths
    syn_data_path = f"../data/synthetic/synthetic_mortality_ground_truth_{beta_form}.csv"
    out_path = f"../data/plotting/final_size_boxplot_data_synthetic_ground_truth_{beta_form}.csv"
    base_res_dir = f"../results/synthetic_ground_truth_{beta_form}"

    # identify zeta grid from the synthetic data columns
    df_head = pd.read_csv(syn_data_path, nrows=1)
    series_cols = [c for c in df_head.columns if c.startswith("y_zeta_")]
    
    # for baseline and the ground-truth variant models
    baseline_model = next(m for m in MODELS if m["name"] == "baseline")
    behavior_model = next(m for m in MODELS if m.get("beta_form") == beta_form)
    models_to_process = [baseline_model, behavior_model]

    rows = []
    series_cols_sorted = sorted(series_cols, key=lambda x: float(x.split('_')[-1]))
    for i, col in enumerate(series_cols_sorted, start=1):
        zeta_val = float(col.split('_')[-1])
        zeta_tag = f"synthetic_zeta{i}"

        for m_dict in models_to_process:
            m_name = m_dict["name"]
            post_path = f"{base_res_dir}/{zeta_tag}/{m_name}_posterior_{RESULT}.csv"

            post_df = pd.read_csv(post_path)
            stats = calculate_extinction_stats(post_df, m_dict)
            stats.update({
                "zeta_true": zeta_val,
                "model_name": m_dict.get("display_name", m_name)
            })
            rows.append(stats)
            
        print(f"Processed {zeta_tag} zeta value = {zeta_val})")

    if rows:
        out_df = pd.DataFrame(rows)
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        out_df.to_csv(out_path, index=False)
        print(f"Saved to {out_path}")

if __name__ == "__main__":
    print("Creating synthetic final size data for boxplots...")
    for variant in ["mixed", "exp", "rational"]:
        generate_final_size_data_for_ground_truth(variant)
