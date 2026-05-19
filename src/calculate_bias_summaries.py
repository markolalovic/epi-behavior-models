#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/calculate_bias_summaries.py

Prints summaries across selected locations for R_0 bias and final size bias.
"""

import numpy as np
import pandas as pd
from config import LOCATIONS, MODELS

def load_medians(path):
    df = pd.read_csv(path)
    df = df[df["location"].isin(LOCATIONS)]
    return df.pivot(index="location", columns="model_name", values="middle").loc[LOCATIONS]

def summarize(x):
    return {
        "n": len(x),
        "n_positive": int(np.sum(x > 0)),
        "median": np.median(x),
        "q25": np.quantile(x, 0.25),
        "q75": np.quantile(x, 0.75),
        "min": np.min(x),
        "max": np.max(x),
    }

def print_summary(label, s, units=""):
    suffix = f" {units}" if units else ""
    print(label)
    print("-" * len(label))
    print(f"positive locations: {s['n_positive']}/{s['n']}")
    print(f"median: {s['median']:.2f}{suffix}")
    print(f"IQR: {s['q25']:.2f} - {s['q75']:.2f}{suffix}")
    print(f"range: {s['min']:.2f} - {s['max']:.2f}{suffix}")
    print()

if __name__ == "__main__":
    baseline = next(m for m in MODELS if not m["is_behavioral"])["display_name"]
    behavior = "Behavioral (Mixed)"

    r0 = load_medians("../data/plotting/r0_boxplot_data.csv")
    final_size = load_medians("../data/plotting/final_size_boxplot_data.csv")

    # baseline underestimates R0, and overestimates final size
    delta_r0 = r0[behavior] - r0[baseline]
    delta_final_pp = 100.0 * (final_size[baseline] - final_size[behavior])

    r0_summary = summarize(delta_r0)
    final_summary = summarize(delta_final_pp)

    print_summary("R0: Behavioral (Mixed) - Baseline", 
                  r0_summary)
    print_summary(
        "Final size: Baseline - Behavioral (Mixed)",
        final_summary,
        units="percentage points"
    )

    max_loc = delta_final_pp.idxmax()
    print("Largest final-size overestimation")
    print(f"location: {max_loc}")
    print(f"difference: {delta_final_pp.loc[max_loc]:.2f} percentage points")

'''
(.venv) ➜  src git:(main) python calculate_bias_summaries.py 
R0: Behavioral (Mixed) - Baseline
---------------------------------
positive locations: 20/20
median: 1.26
IQR: 0.91 - 1.97
range: 0.15 - 2.67

Final size: Baseline - Behavioral (Mixed)
-----------------------------------------
positive locations: 20/20
median: 2.96 percentage points
IQR: 1.58 - 6.22 percentage points
range: 0.64 - 7.87 percentage points

Largest final-size overestimation
location: WI
difference: 7.87 percentage points
'''