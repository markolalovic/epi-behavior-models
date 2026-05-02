#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/utils.py
Some math and data helpers.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tabulate import tabulate
from config import SYNTHETIC_CONFIG, RESULT

def trailing_ma7(x):
    # x: 1D np.array of length T
    k = 7
    kernel = np.ones(k, dtype=float) / k
    if len(x) < k:
        # fill with 0
        return np.zeros_like(x, dtype=float)
    valid = np.convolve(x, kernel, mode="valid")  # length T - k + 1
    return np.concatenate([np.zeros(k - 1, dtype=float), valid])

def weighted_quantile(values, weights, quantile=0.5):
    """
    Compute a weighted quantile:
     * sort the values and accumulate their normalized weights
     * compute the cumulative sum (CDF)
     * see what fraction of the total weight is below each value
     * then interpolate the quantile (0.5 for the median) 
     * along this cumulative to get the corresponding value
    """
    values = np.asarray(values)
    weights = np.asarray(weights)

    order = np.argsort(values)
    values, weights = values[order], weights[order] / np.sum(weights)

    cdf = np.cumsum(weights)
    cdf /= cdf[-1]

    return np.interp(quantile, cdf, values)

def variant_name(beta_form):
    mapping = {
        "exp": "behavior_exp",
        "mixed": "behavior_mixed",
        "rational": "behavior_rational",
    }
    return mapping[beta_form]

def ground_truth_values(beta_form=None, baseline=False, gt_location=None):
    """Extracts empirically plausible parameter values for given model 
    and specified ground truth location in config.SYNTHETIC_CONFIG."""
    if baseline:
        model_name = "baseline"
    else:
        model_name = variant_name(beta_form)
    
    if not gt_location:
        gt_location = SYNTHETIC_CONFIG["ground_truth_location"]
    posterior_path = f"../results/{gt_location}/{model_name}_posterior_{RESULT}.csv"
    df = pd.read_csv(posterior_path)

    weights = df["weight"].to_numpy()
    weights = weights / np.sum(weights)

    median_values = {
        "theta_pi0": weighted_quantile(df["theta_pi0"], weights),
        "R0": weighted_quantile(df["R0"], weights),
        "delta": weighted_quantile(df["delta"], weights)
    }
    if not baseline:
        median_values["zeta"] = weighted_quantile(df["zeta"].to_numpy(), weights)
    return median_values

def extract_parameters(location):
    """Extracts empirically plausible parameter values that 
    can be used as ground-truth values.
    """
    def print_params(name, params):
        print(f"{name} = {{")
        for k, v in params.items():
            print(f"    '{k}': {v:.6f},")
        print("}")
    
    print(f"Extracted weighted medians for ground truth location {location}:")
    base = ground_truth_values(baseline=True, gt_location=location)
    print_params("EXAMPLE_BASELINE", base)

    mixed = ground_truth_values(beta_form="mixed", baseline=False, gt_location=location)
    print_params("EXAMPLE_BEHAVIOR", mixed)    

def print_nice_table(df):
    """Prints a nice Markdown formated table with the minimum values in bold."""
    df_fmt = df.copy().astype(object)
    numeric_cols = [c for c in df.columns if c != "Location"]
        
    for i, row in df.iterrows():
        min_val = row[numeric_cols].min()
        for col in numeric_cols:
            val = row[col]
            if val == min_val:
                df_fmt.at[i, col] = f"**{val:.4f}**"
            else:
                df_fmt.at[i, col] = f"{val:.4f}"
    print(
        tabulate(df_fmt, 
                 headers='keys', 
                 tablefmt="github", 
                 stralign="center", 
                 numalign="right"
                 )
    )
