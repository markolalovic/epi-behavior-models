#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/inference.py
ABC-SMC inference pipeline:
- distance function
- run_abc function
"""

import os
import numpy as np
import pickle as pkl
from datetime import timedelta
from pyabc import ABCSMC, MultivariateNormalTransition, QuantileEpsilon
from config import ABCSMC_CONFIG, ABC_RUN_CONFIG, RESULT

def nssr_distance(simulation_output, observation_dict):
    r"""
    Normalized sum of squared error:
      sum((y_sim - y_obs)^2) / sum(y_obs^2)
    """
    sim_data = simulation_output.get("data")
    obs_data = observation_dict.get("data")

    if sim_data is None or np.isinf(sim_data).any():
        return np.inf

    obs_data = np.asarray(obs_data)
    sim_data = np.asarray(sim_data)
    
    numerator = np.sum(np.square(sim_data - obs_data))
    denominator = np.sum(np.square(obs_data))

    if denominator == 0:
        return numerator
    
    return numerator / denominator

def run_abc(fixed_params, obs_deaths, base_prior, location, model_dict):
    """
    Runs the ABC-SMC inference for a specific location and model
    """
    model_name = model_dict['name']
    print(f"\n--- Starting ABC-SMC for: {model_name} ({location}) ---")

    from models import run_seird_model
    
    def model_wrapper(params):
        return run_seird_model(
            params=params,
            fixed_params=fixed_params,
            obs_length=len(obs_deaths),
            model_dict=model_dict
        )

    params_prior = base_prior.copy()
    if not model_dict['is_behavioral']:
        params_prior.pop("zeta")

    # set ABC-SMC using settings from config.py
    transition_kernel = MultivariateNormalTransition()

    # aggressive enough to focus the posterior
    eps_sched = QuantileEpsilon(alpha=ABCSMC_CONFIG["quantile_alpha"])

    abc = ABCSMC(
        models=model_wrapper,
        parameter_priors=params_prior,
        distance_function=nssr_distance,
        population_size=ABCSMC_CONFIG["population_size"],
        transitions=transition_kernel,
        eps=eps_sched
    )

    # set up persistent database path
    db_dir = "../data/local-db/"
    os.makedirs(db_dir, exist_ok=True)

    # NOTE: avoiding path separators for DB filename
    # No nested directories
    # use RESULT flag in the filename to keep runs separate
    safe_loc = str(location).replace(os.sep, "_").replace("/", "_")    
    db_path = os.path.join(db_dir, f"{safe_loc}_{model_name}_{RESULT}.db")
    if os.path.exists(db_path):
        os.remove(db_path)
    
    # --- run the inference ---
    abc.new("sqlite:///" + db_path, {"data": obs_deaths})

    max_walltime_s = ABC_RUN_CONFIG["max_walltime_s"]
    history = abc.run(
        max_nr_populations=ABC_RUN_CONFIG["max_nr_populations"],
        max_walltime=timedelta(seconds=max_walltime_s),
        max_total_nr_simulations=ABC_RUN_CONFIG["max_total_nr_simulations"],
        min_eps_diff=ABC_RUN_CONFIG["min_eps_diff"]
    )

    # --- extract the results ---
    results_dir = f"../results/{location}/"
    os.makedirs(results_dir, exist_ok=True)
    
    history_path = os.path.join(results_dir, f"{model_name}_history_{RESULT}.pkl")
    with open(history_path, 'wb') as file:
        pkl.dump(history, file)
    print(f"Saved full history to: {history_path}")

    df, weights = history.get_distribution(m=0, t=history.max_t)
    df["weight"] = weights
    csv_path = os.path.join(results_dir, f"{model_name}_posterior_{RESULT}.csv")
    df.to_csv(csv_path, index=False)
    print(f"Saved final posterior samples and weights to: {csv_path}")
    
    return history
