#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/create_fits_data.py
Generates the data needed for plotting the model fits figures:
- Main Figure 1, and Supplement Figures S4, S5

For each location and model:
- Loads the samples and weights of estimated posterior distribution
- Samples 1000 parameter sets from the posterior using weights
- Runs a simulation for each sampled parameter set
- Calculates the median and 90% predictive interval, i.e. 0.05, 0.5, 0.95 quantiles
- Aggregates all results and save them to a single CSV file for figures: `data/plotting/fits_data.csv`

"""

import os
import pandas as pd
import numpy as np
from config import MODELS, LOCATIONS, POPULATION_SIZE, FIXED_PARAMS, RESULT
from models import run_seird_model

def generate_predictive_trajectories(posterior_df, model_dict, location, obs_length, n_samples=1000):
    """
    Samples from estimated posterior to generate ensamble of simulated trajectories.
    """
    print(f"Sampling {n_samples} trajectories...")
    # sample parameter sets from the posterior using weights
    sampled_params = posterior_df.sample(
        n=n_samples,
        weights='weight',
        replace=True,
        random_state=None
    ).to_dict('records')

    # prepare location-specific fixed params
    fixed_params = FIXED_PARAMS.copy()
    fixed_params["N"] = POPULATION_SIZE[location]

    # simulate for each draw
    trajectories = []
    for params in sampled_params:
        sim_out = run_seird_model(
            params=params,
            fixed_params=fixed_params,
            obs_length=obs_length,
            model_dict=model_dict
        )
        trajectories.append(sim_out['data'])

    return np.array(trajectories)


def generate_fits_data():
    """
    Generate and save the data for the model fits figure.
    """
    print("Generating data for model fits figure")

    models_to_process = MODELS

    # load the observed data
    data_path = "../data/processed/smoothed_mortality.csv"
    df_all = pd.read_csv(data_path, parse_dates=["date"])
    df_all = df_all.sort_values("date").reset_index(drop=True)
    T = len(df_all)
    print(f"Loaded observed data from {data_path}, of length {T} days.")

    all_results = []
    for location in LOCATIONS:
        print(f"\nProcessing location: {location}")

        # data series for this location is already smoothed
        obs_series = df_all[location].to_numpy(dtype=float)
        obs_length = len(obs_series)

        for model_dict in models_to_process:
            model_name = model_dict['name']
            posterior_path = f"../results/{location}/{model_name}_posterior_{RESULT}.csv"

            print(f"Processing model: {model_name}")
            posterior_df = pd.read_csv(posterior_path)

            # sample and simulate to get 1000 trajectories
            trajectories = generate_predictive_trajectories(
                posterior_df=posterior_df,
                model_dict=model_dict,
                location=location,
                obs_length=obs_length,
                n_samples=1000
            )

            # calculate the quantiles
            q_low, q_med, q_hi = np.quantile(trajectories, [0.05, 0.5, 0.95], axis=0)

            # 5) Assemble tidy block for this (location, model)
            block = pd.DataFrame({
                'date': df_all['date'].values,                  # shared date axis
                'location': location,
                'model_name': model_dict.get('display_name', model_name),
                'lower_90': q_low,
                'median': q_med,
                'upper_90': q_hi,
                'obs': obs_series
            })

            all_results.append(block)

    final_df = pd.concat(all_results, ignore_index=True)

    # save the aggregated data to a csv file
    out_dir = "../data/plotting/"
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "fits_data.csv")
    final_df.to_csv(out_path, index=False)
    print(f"\nSaved model fits data to: {out_path}")


if __name__ == "__main__":
    generate_fits_data()
