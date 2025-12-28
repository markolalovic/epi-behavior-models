#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""create_r0_data.py
Calculates the statistics for R_0 boxplots Panel A figures:
- Main Figure 2, and Supplement Figures S6, S7

Outputs:
- `data/plotting/r0_boxplot_data.csv`
"""

import os
import pandas as pd
import numpy as np
from statsmodels.stats.weightstats import DescrStatsW
from config import MODELS, LOCATIONS, RESULT

def weighted_quartiles(x, w):
    """
    Calculates weighted quartiles using statsmodels DescrStatsW
    Returns: q1, median, q3
    """
    w = np.asarray(w, dtype=float)
    x = np.asarray(x, dtype=float)
    # normalize the weights, have to be strictly positive
    w = w / w.sum() if w.sum() > 0 else np.full_like(w, 1.0 / len(w))
    ds = DescrStatsW(x, weights=w, ddof=0)
    q1, med, q3 = ds.quantile([0.25, 0.5, 0.75], return_pandas=False)
    return (float(q1), float(med), float(q3))

def boxplot_stats_from_weighted_samples(x, w):
    """
    Calculates boxplot stats from weighted samples:
      - quartiles are weighted
      - whiskers use 1.5 * IQR rule applied to x, and min / max within bounds
    Returns dict with keys: ymin, lower, middle, upper, ymax
    """
    q1, med, q3 = weighted_quartiles(x, w)
    iqr = q3 - q1
    lower_lim = q1 - 1.5 * iqr
    upper_lim = q3 + 1.5 * iqr

    # whiskers: closest points within [lower_lim, upper_lim]
    x_in = x[(x >= lower_lim) & (x <= upper_lim)]
    ymin = np.min(x_in)
    ymax = np.max(x_in)
    return dict(ymin=ymin, lower=q1, middle=med, upper=q3, ymax=ymax)


def generate_r0_boxplot_data():
    """
    Generates and saves the R0 boxplot data.
    """
    print("Generating data for R0 boxplots...")

    models_to_process = MODELS
    rows = []
    for location in LOCATIONS:
        for model in models_to_process:
            model_name = model['name']
            display_name = model.get('display_name', model_name)
            posterior_path = f"../results/{location}/{model_name}_posterior_{RESULT}.csv"

            df = pd.read_csv(posterior_path)

            r0_vals = df['R0'].to_numpy(dtype=float)
            w = df['weight'].to_numpy(dtype=float)
            
            # normalize weights, they have to be nonnegative
            w = w / w.sum() if w.sum() > 0 else np.full_like(w, 1.0 / len(w))

            stats = boxplot_stats_from_weighted_samples(r0_vals, w)
            stats['location'] = location
            stats['model_name'] = display_name
            rows.append(stats)

    out_df = pd.DataFrame(rows)
    out_dir = "../data/plotting/"
    out_path = os.path.join(out_dir, "r0_boxplot_data.csv")
    out_df.to_csv(out_path, index=False)
    print(f"\nSuccessfully saved R0 boxplot data to: {out_path}")

if __name__ == "__main__":
    generate_r0_boxplot_data()
