#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/create_synthetic_r0_data.py
Calculates data for synthetic experiments Panel B figures:
- Weighted boxplot statistics for R_0
- Weighted quartiles and 1.5 * IQR rule for whiskers

For each ground truth beta_form in ["mixed", "exp", "rational"]:
  For each synthetic_zeta{idx} and for each model {baseline, behavioral with that beta_form}:
  - Loads synthetic smoothed mortality:
    `data/synthetic/synthetic_mortality_ground_truth_{beta_form}.csv`
  - Loads results (columns: R0, and weight) 
    `results/synthetic_ground_truth_{beta_form}/synthetic_zeta{idx}/{model}_posterior_{RESULT}.csv`
  - Computes weighted quartiles (Q1, median, Q3) and whiskers using 1.5*IQR rule.

Outputs:
- data/plotting/r0_boxplot_data_synthetic_ground_truth_mixed.csv
- data/plotting/r0_boxplot_data_synthetic_ground_truth_exp.csv
- data/plotting/r0_boxplot_data_synthetic_ground_truth_rational.csv

"""

import numpy as np
import pandas as pd
from config import MODELS, RESULT, SYNTHETIC_CONFIG
from utils import weighted_quantile, variant_name

def get_weighted_box_stats(x, w):
    """Calculates weighted boxplot statistics: ymin, lower, middle, upper, ymax."""
    # weighted quartiles
    q1, med, q3 = weighted_quantile(x, w, [0.25, 0.5, 0.75])
    
    # whiskers using 1.5 * IQR rule
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    
    # whiskers are the min / max values within the bounds
    xin = x[(x >= lower_bound) & (x <= upper_bound)]
    if xin.size == 0:
        return {"ymin": x.min(), "lower": q1, "middle": med, "upper": q3, "ymax": x.max()}
    
    return {
        "ymin": float(np.min(xin)),
        "lower": float(q1),
        "middle": float(med),
        "upper": float(q3),
        "ymax": float(np.max(xin))
    }

def generate_r0_data_for_ground_truth(beta_form):
    """G R0 posteriors for a specific behavioral ground truth."""
    print(f"\nCreating synthetic R0 data | Ground Truth: {beta_form}")

    # data paths
    syn_data_path = f"../data/synthetic/synthetic_mortality_ground_truth_{beta_form}.csv"
    out_path = f"../data/plotting/r0_boxplot_data_synthetic_ground_truth_{beta_form}.csv"
    base_res_dir = f"../results/synthetic_ground_truth_{beta_form}"

    # load synthetic data header to get zeta values
    df_wide = pd.read_csv(syn_data_path, nrows=1)
    series_cols = [c for c in df_wide.columns if c.startswith("y_zeta_")]
    
    # process all models: baseline and the variant used for truth
    baseline_model = next(m for m in MODELS if m["name"] == "baseline")
    behavior_model = next(m for m in MODELS if m.get("beta_form") == beta_form)
    models_to_process = [baseline_model, behavior_model]

    rows = []
    # for all 6 zeta values in a grid
    series_cols_sorted = sorted(series_cols, key=lambda x: float(x.split('_')[-1]))
    for i, col in enumerate(series_cols_sorted, start=1):
        zeta_val = float(col.split('_')[-1])
        zeta_tag = f"synthetic_zeta{i}"

        for m_dict in models_to_process:
            m_name = m_dict["name"]
            post_path = f"{base_res_dir}/{zeta_tag}/{m_name}_posterior_{RESULT}.csv"

            df = pd.read_csv(post_path)
            x = df["R0"].to_numpy(dtype=float)
            w = df["weight"].to_numpy(dtype=float)

            stats = get_weighted_box_stats(x, w)
            stats.update({
                "zeta_true": zeta_val,
                "model_name": m_dict.get("display_name", m_name)
            })
            rows.append(stats)

    out_df = pd.DataFrame(rows)
    out_df.to_csv(out_path, index=False)
    print(f"Saved to {out_path}")

if __name__ == "__main__":
    print("Creating synthetic R_0 data for boxplots...")
    for variant in ["mixed", "exp", "rational"]:
        generate_r0_data_for_ground_truth(variant)
