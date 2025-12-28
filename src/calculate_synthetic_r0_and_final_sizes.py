#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""src/calculate_synthetic_r0_and_final_sizes.py
Calculates reference ground-truth values for R_0 and final epidemic sizes.

Ground-truth R_0 values:
  - are weighted medians from posterior for MA 
  - they depend behavior variant beta_form
  - so 3 values in total

Ground-truth final epidemic sizes:
  - are estimated from running the models 
  - as specified in compute_final_size_until_extinction()
  - they depend on behavior variant beta_form,
  - and on zeta values as specified in SYNTHETIC_CONFIG["zeta_grid"]
  - so 3 * 6 = 18 values in total

Outputs:
- data/plotting/synthetic_true_r0.csv
- data/plotting/synthetic_true_final_sizes.csv
"""

import os
import pandas as pd

from config import (
    MODELS, POPULATION_SIZE, FIXED_PARAMS, SYNTHETIC_CONFIG
)

from models import compute_final_size_until_extinction
from utils import ground_truth_values

# config
zeta_grid = SYNTHETIC_CONFIG["zeta_grid"]
gt_location = SYNTHETIC_CONFIG["ground_truth_location"]
out_dir = "../data/plotting/"

if __name__ == "__main__":
    r0_data = {}
    final_sizes = {"zeta_true": zeta_grid}
    
    fixed_params = FIXED_PARAMS.copy()
    fixed_params["N"] = POPULATION_SIZE[gt_location]

    print(f"Calculating ground truth values for {gt_location}")

    for variant in ["mixed", "exp", "rational"]:
        model_name = f"behavior_{variant}"
        model_dict = next(m for m in MODELS if m["name"] == model_name)
        
        # medians from gt_location posterior using helper from utils.py
        base_params = ground_truth_values(variant)
        r0_data[f"true_r0_{variant}"] = [base_params["R0"]]
        
        # calculate final sizes over for each zeta
        variant_final_sizes = []
        print(f"Calculating {variant} final sizes...")
        for z in zeta_grid:
            params_z = dict(base_params)
            params_z["zeta"] = z
            
            # run simulation until extinction
            _, _, fs = compute_final_size_until_extinction(
                params=params_z,
                fixed_params=fixed_params,
                model_dict=model_dict,
                total_days=1825   # run for ~ 5 years just in case
            )
            variant_final_sizes.append(fs)
        
        final_sizes[f"true_final_size_{variant}"] = variant_final_sizes

    r0_path = os.path.join(out_dir, "synthetic_true_r0.csv")
    pd.DataFrame(r0_data).to_csv(r0_path, index=False) 
    print(f"Saved true R_0 values to: {r0_path}")
    
    fs_path = os.path.join(out_dir, "synthetic_true_final_sizes.csv")
    pd.DataFrame(final_sizes).to_csv(fs_path, index=False)
    print(f"Saved true final sizes to: {fs_path}")
