#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" run_inference.py
Runs the ABC-SMC inference for all specified models and locations.
"""

import time
import numpy as np
import pandas as pd

from config import LOCATIONS, MODELS, POPULATION_SIZE
from config import FIXED_PARAMS, PRIORS

from inference import run_abc

if __name__ == "__main__":
    print(f"--- selected {len(LOCATIONS)} states: ----")
    print(LOCATIONS)
    print("-------------------------")

    start_time = time.time()
    print("--- Starting all ABC-SMC fitting runs ---")

    data_path = "../data/smoothed_mortality.csv"
    df_all = pd.read_csv(data_path, parse_dates=["date"])
    df_all = df_all.sort_values("date").reset_index(drop=True)

    for location in LOCATIONS:
        print(f"\n==================================================")
        print(f"  Processing location: {location}")
        print(f"==================================================")

        # extract the smoothed observed series for this location
        obs_deaths = df_all[location].to_numpy(dtype=float)
        T = len(obs_deaths)
        print(f"Loaded observed data from {data_path}. Length: {T} days.")

        # prepare location-specific fixed parameters
        current_fixed_params = FIXED_PARAMS.copy()
        current_fixed_params["N"] = POPULATION_SIZE[location]

        # loop over models
        for model_dict in MODELS:
            run_abc(
                fixed_params=current_fixed_params,
                obs_deaths=obs_deaths,
                base_prior=PRIORS,
                location=location,
                model_dict=model_dict
            )

    elapsed = (time.time() - start_time) / 60.0
    print(f"\n==================================================")
    print(f"All runs completed in {elapsed:.2f} minutes.")
    print(f"==================================================")
