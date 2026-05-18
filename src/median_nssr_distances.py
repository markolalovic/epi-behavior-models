#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/median_nssr_distances.py
Computes weighted median posterior NSSE distances for all models and locations.

Outputs:
- `results/model_comparison/median_nssr_distances.csv`
"""

import os
import numpy as np
import pandas as pd
from tabulate import tabulate

from config import (
    MODELS, POPULATION_SIZE, FIXED_PARAMS, 
    RESULT, LOCATIONS, LOCATION_NAME, LOCATIONS_ALL
)

from models import run_seird_model
from inference import nssr_distance
from utils import weighted_quantile, print_nice_table

if __name__ == "__main__":
    obs_path = "../data/processed/smoothed_mortality.csv"
    obs_df = pd.read_csv(obs_path, parse_dates=["date"]).sort_values("date").reset_index(drop=True)

    print(f"Calculating weighted median posterior distances...")

    rows = []
    for loc in LOCATIONS_ALL:
        print(f"Processing {LOCATION_NAME[loc]}...")
        obs = obs_df[loc].to_numpy(dtype=float)
        obs_dict = {"data": obs}
        
        row = {"Location": LOCATION_NAME[loc], "Abbr": loc}

        # for each model: baseline and 3 behavior variants:
        for model in MODELS:
            model_name = model['name']
            post_path = f"../results/{loc}/{model_name}_posterior_{RESULT}.csv"

            if not os.path.exists(post_path):
                row[model['display_name']] = np.nan
                continue

            df = pd.read_csv(post_path)
            weights = df["weight"].to_numpy()
            
            distances = []
            fixed_params = FIXED_PARAMS.copy()
            fixed_params["N"] = POPULATION_SIZE[loc]

            # simulate for each particle to get the distance
            for _, params in df.iterrows():
                params_dict = {k: float(params[k]) for k in df.columns if k != "weight"}
                
                sim_out = run_seird_model(
                    params=params_dict,
                    fixed_params=fixed_params,
                    obs_length=len(obs),
                    model_dict=model
                ) 
                distances.append(nssr_distance(sim_out, obs_dict))
            # calc weighted median distance
            row[model['display_name']] = weighted_quantile(distances, weights, 0.5)
        
        # save the results for each model
        rows.append(row)

    # show results
    df_out = pd.DataFrame(rows)
    print("Weighted Median NSSR")
    print_nice_table(df_out.drop(columns=["Abbr"]))

    # export results
    out_dir = "../results/model_comparison"
    out_path = os.path.join(out_dir, f"median_nssr_distances_{RESULT}.csv")
    df_out.to_csv(out_path, index=False)
    print(f"\nSaved median NSSR distances to: {out_path}")
