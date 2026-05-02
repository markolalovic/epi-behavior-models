#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/create_final_epi_size_data.py
Calculates the statistics for final epidemic size boxplots Panel B figures:
- Main Figure 2, and Supplement Figures S6, S7

The final epidemic size is calculated from estimated posterior distributions as:

    final_size = 1 - S(t_end) / N

where t_end is the first day such that I(t_end) < 1.0, obtained by
`compute_final_size_until_extinction()` from `src/models.py`.

Outputs:
- `data/plotting/final_size_boxplot_data.csv`
"""

import os
import pandas as pd
import numpy as np
from config import MODELS, LOCATIONS, POPULATION_SIZE, FIXED_PARAMS, RESULT
from models import compute_final_size_until_extinction

def calculate_final_size_distribution(posterior_df, model_dict, location,
    n_samples = 1000, total_days = 1825):
    """
    Calculates posterior predictive distribution of the **extinction-based**
    final epidemic size:
     total_days = 1825 ~ 5 years
# integrate long enough to let I(t) die out for *real* data too
        final_size = 1 - S(t_end) / N

    by:
      1. drawing n_samples parameter sets from the weighted posterior,
      2. running the model to (up to) total_days,
      3. extracting S(t_end) from compute_final_size_until_extinction(...).
    """
    print(f"  Generating {n_samples} extinction-based final sizes...")

    # 1) sample params (weighted)
    sampled_params = posterior_df.sample(
        n=n_samples, weights="weight", replace=True
    ).to_dict("records")

    # 2) fixed params for this location
    fixed_params = FIXED_PARAMS.copy()
    fixed_params["N"] = POPULATION_SIZE[location]
    N = float(fixed_params["N"])

    final_sizes = []
    for params in sampled_params:
        # weight column is already dropped by to_dict
        # → but if the CSV had extra columns, be defensive:
        if "weight" in params:
            params = {k: v for k, v in params.items() if k != "weight"}

        # run to extinction
        t_end, S_t_end, fin = compute_final_size_until_extinction(
            params=params,
            fixed_params=fixed_params,
            model_dict=model_dict,
            total_days=total_days,
        )
        # `fin` should already be 1 - S(t_end)/N, but keep a safety check:
        if not np.isfinite(fin):
            # fall back to computing from S_t_end
            fin = 1.0 - (S_t_end / N)

        final_sizes.append(fin)

    return np.asarray(final_sizes, dtype=float)


def get_boxplot_stats_from_sample(sample: np.ndarray) -> dict:
    """
    Calculates boxplot statistics from an unweighted sample.
    Returns a dict with keys compatible with ggplot's stat="identity".
    """
    sample = np.asarray(sample, dtype=float)
    if sample.size == 0:
        return {
            "ymin": np.nan,
            "lower": np.nan,
            "middle": np.nan,
            "upper": np.nan,
            "ymax": np.nan,
        }

    q1, median, q3 = np.quantile(sample, [0.25, 0.5, 0.75])
    iqr = q3 - q1
    lower_whisker_limit = q1 - 1.5 * iqr
    upper_whisker_limit = q3 + 1.5 * iqr

    # clip whiskers to points inside the limits
    whisker_low = sample[sample >= lower_whisker_limit].min()
    whisker_high = sample[sample <= upper_whisker_limit].max()

    return {
        "ymin": float(whisker_low),
        "lower": float(q1),
        "middle": float(median),
        "upper": float(q3),
        "ymax": float(whisker_high),
    }

def generate_final_size_data():
    """
    Main entry point.
    Loads all posteriors, recomputes final epidemic size with the *new*
    definition, and saves one CSV for plotting.
    """
    print("--- Generating data for Final Epidemic Size boxplots (extinction-based) ---")

    models_to_process = MODELS

    # 0) Load observed, wide
    data_path = "../data/processed/smoothed_mortality.csv"
    df_all = pd.read_csv(data_path, parse_dates=["date"]).sort_values("date").reset_index(drop=True)
    T = len(df_all)
    print(f"Loaded observed data from {data_path}. Observation length: {T} days.")

    # Sanity: check locations exist
    missing = [loc for loc in LOCATIONS if loc not in df_all.columns]
    if missing:
      raise RuntimeError(f"Missing columns in smoothed_mortality.csv: {missing}")

    rows = []

    for location in LOCATIONS:
        print(f"\nLocation: {location}")

        for model_dict in models_to_process:
            model_name = model_dict["name"]
            disp_name = model_dict.get("display_name", model_name)
            post_path = f"../results/{location}/{model_name}_posterior_{RESULT}.csv"

            if not os.path.exists(post_path):
                print(f"  Warning: missing posterior for {model_name} @ {location}: {post_path}")
                continue

            posterior_df = pd.read_csv(post_path)
            if "weight" not in posterior_df.columns:
                print(f"  Warning: 'weight' missing in {post_path}; skipping.")
                continue

            # compute extinction-based sample
            sample = calculate_final_size_distribution(
                posterior_df=posterior_df,
                model_dict=model_dict,
                location=location,
                n_samples=1000,
                total_days=1825,
            )

            stats = get_boxplot_stats_from_sample(sample)
            stats["location"] = location
            stats["model_name"] = disp_name
            rows.append(stats)

    out_df = pd.DataFrame(rows)
    out_dir = "../data/plotting/"
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "final_size_boxplot_data.csv")
    out_df.to_csv(out_path, index=False)

    print(f"\nSaved extinction-based final size data to: {out_path}")


if __name__ == "__main__":
    generate_final_size_data()
