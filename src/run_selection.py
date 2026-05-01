#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/run_model_selection.py
Performs ABC-SMC model selection between the Baseline model 
and each of the three behavioral variants: mixed, exp, and rational.

Outputs:
- Summary results for each variant across 30 locations: 
  * `results/model_selection/baseline_vs_{variant_name}_summary_{RESULT}.csv`

- Convergence plots of posterior model probabilities P(M | D, t) over generations: 
  * `results/model_selection/{location}_{behavior_name}_trajectory_{RESULT}.png`

- Raw pyABC history objects: particles, distances and weights for each generation: 
  * `data/local-db/{location}_selection_{behavior_name}_{RESULT}.db`

"""

import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tabulate import tabulate

from pyabc import ABCSMC, MultivariateNormalTransition, QuantileEpsilon
from config import (
    LOCATIONS, MODELS, POPULATION_SIZE, FIXED_PARAMS, 
    RESULT, LOCATION_NAME, PRIORS, MODEL_SELECTION_CONFIG
)

from models import run_seird_model
from inference import nssr_distance

def prior_for_model(model_dict):
    """Behavioral variants use full prior, baseline has no zeta."""
    if model_dict["is_behavioral"]:
        return PRIORS
    else:
        base_prior = PRIORS.copy()
        if "zeta" in base_prior:
            base_prior.pop("zeta")
        return base_prior

def extract_model_prob_trajectory(history, model_names):
    """Returns dataframe with 
      `t, baseline, behavior(variant)`
    where each column is posterior model probability at generation t.
    """
    rows = []
    for t in range(history.max_t + 1):
        s = history.get_model_probabilities(t=t)  
        row = {"t": t}
        for i, name in enumerate(model_names):
            row[name] = float(s.loc[i, "p"]) if i in s.index else 0.0
        rows.append(row)
    return pd.DataFrame(rows)

def find_stable_generation(df, behavior_name, baseline_name, delta, k_consec):
    """
    Picks the first `t` where: 
    - proabilities change < delta 
    - for k_conseq generations 
    - and P(Baseline) > 0
    """
    T = int(df["t"].max())
    df = df.sort_values("t").reset_index(drop=True)
    df["dB"] = df[behavior_name].diff().abs()
    df["dBase"] = df[baseline_name].diff().abs()

    for t in range(1, T + 1):
        t0 = max(1, t - k_consec + 1)
        window = df[(df["t"] >= t0) & (df["t"] <= t)]
        if len(window) == k_consec and (window["dB"].max() < delta) and (window["dBase"].max() < delta):
            if float(df.loc[df["t"] == t, baseline_name].iloc[0]) > 0.0:
                return t
    
    # else, fallback to last known generation where baseline was not extinct yet
    nonzero = df.loc[df[baseline_name] > 0.0, "t"]
    if len(nonzero) > 0:
        return int(nonzero.iloc[-1])
    return T

def equal_two_model_bayes(p_behavior, p_baseline):
    r"""Calculates BF_eq under equal model priors:
      BF = P(behavior | D) / P(baseline | D)
      P_equal(behavior | D) = BF / (1 + BF) 
    """
    if p_baseline <= 1e-12:
        return np.inf, 1.0
    bf = p_behavior / p_baseline
    p_equal = bf / (1.0 + bf)
    return bf, p_equal

def interpret_bf(bf):
    """Interpret Bayes factor behavior vs baseline. """
    if not np.isfinite(bf): return "very strong (Behavior)"
    if bf < 1: return "favors Baseline"
    if 1 <= bf < 3: return "very weak"
    if 3 <= bf < 20: return "positive"
    if 20 <= bf < 150: return "strong"
    return "very strong"

def run_selection_pair(location, behavior_model_cfg, 
                       df_all, pop_size, max_pops, alpha_q, stab_delta, stab_k):
    """Runs model selection for a single location from baseline and behavior(variant)."""

    behavior_name = behavior_model_cfg["name"]
    baseline_cfg = next(m for m in MODELS if m["name"] == "baseline")
    
    print(f"\nRunning model selection: Baseline vs {behavior_name} for location: {location}")    
    
    model_configs = [baseline_cfg, behavior_model_cfg]
    model_names = [m["name"] for m in model_configs]
    
    # observed data
    obs_deaths = df_all[location].to_numpy()
    obs_length = len(obs_deaths)

    # location specific fixed parameters
    fixed_params_loc = FIXED_PARAMS.copy()
    fixed_params_loc["N"] = POPULATION_SIZE[location]

    # model functions and priors
    model_fns = []
    priors = []
    transitions = []
    for m in model_configs:
        def model_fn(params, model_config=m):
            return run_seird_model(
                params=params, 
                fixed_params=fixed_params_loc, 
                obs_length=obs_length, 
                model_dict=model_config
            )
        model_fns.append(model_fn)
        priors.append(prior_for_model(m))
        transitions.append(MultivariateNormalTransition())

    # run ABC-SMC
    abc = ABCSMC(
        models=model_fns, 
        parameter_priors=priors, 
        distance_function=nssr_distance,
        population_size=pop_size, 
        transitions=transitions, 
        eps=QuantileEpsilon(alpha=alpha_q)
    )

    # set DB path for the current pairwise run
    db_path = f"../data/local-db/{location}_selection_{behavior_name}_{RESULT}.db"
    if os.path.exists(db_path):
        os.remove(db_path)

    abc.new("sqlite:///" + db_path, {"data": obs_deaths})
    history = abc.run(max_nr_populations=max_pops)

    # extract and save results
    out_dir = "../results/model_selection/"
    os.makedirs(out_dir, exist_ok=True)

    df_traj = extract_model_prob_trajectory(history, model_names)
    
    # stabilization hyperparameters
    t_star = find_stable_generation(
        df_traj, 
        behavior_name, 
        baseline_cfg["name"], 
        delta=stab_delta, k_consec=stab_k
    )
    
    row_star = df_traj.loc[df_traj["t"] == t_star].iloc[0]
    bf_eq, _ = equal_two_model_bayes(row_star[behavior_name], row_star[baseline_cfg["name"]])

    # plot model selection convergence
    plt.figure(figsize=(7, 4))
    plt.plot(df_traj["t"], df_traj[behavior_name], marker="o", label=behavior_model_cfg["display_name"])
    plt.plot(df_traj["t"], df_traj[baseline_cfg["name"]], marker="o", label="Baseline")
    plt.axvline(t_star, color="red", linestyle="--", label=f"g* = {t_star}")
    plt.xlabel("Generation t")
    plt.ylabel("P(model | D, t)")
    plt.ylim(0, 1.05)
    plt.legend()
    plt.title(f"{LOCATION_NAME[location]}: model selection convergence")
    plt.savefig(f"{out_dir}{location}_{behavior_name}_trajectory_{RESULT}.png")
    plt.close()

    return {
        "location": location, 
        "behavior_variant": behavior_name, 
        "t_star": t_star,
        "P_baseline": row_star[baseline_cfg["name"]], 
        "P_behavior": row_star[behavior_name],
        "BF_eq": bf_eq, 
        "Interpretation": interpret_bf(bf_eq)
    }

if __name__ == "__main__":
    t0 = time.time()

    # hyperparameters for model selection
    pop_size = MODEL_SELECTION_CONFIG["population_size"]
    max_pops = MODEL_SELECTION_CONFIG["max_nr_populations"]
    alpha_q  = MODEL_SELECTION_CONFIG["quantile_alpha"]
    stab_delta = MODEL_SELECTION_CONFIG["stab_delta"]
    stab_k     = MODEL_SELECTION_CONFIG["stab_k"]

    # observed smoothed mortality data for all locations
    data_path = "../data/processed/smoothed_mortality.csv"
    df_all = pd.read_csv(data_path, parse_dates=["date"]) \
               .sort_values("date").reset_index(drop=True)

    # for all behavioral variants vs baseline - i.e. three pairwise comparisons!
    behavioral_variants = [m for m in MODELS if m["is_behavioral"]]
    
    for variant in behavioral_variants:
        # variant name for filename, only: 'mixed', 'exp', 'rational'
        variant_name = variant['name'].replace("behavior_", "")
        
        print(f"\n{'='*60}")
        print(f"Model selection: Baseline vs {variant['display_name']}")
        print(f"{'='*60}")

        variant_results = []
        for loc in LOCATIONS:
            # run the selection for the pair for this location
            res = run_selection_pair(
                location=loc, 
                behavior_model_cfg=variant, 
                df_all=df_all,
                pop_size=pop_size, 
                max_pops=max_pops, 
                alpha_q=alpha_q, 
                stab_delta=stab_delta, 
                stab_k=stab_k
            )
            variant_results.append(res)
        
        df_variant = pd.DataFrame(variant_results)
        print(f"\nSummary for Baseline vs {variant['display_name']}:")
        print(tabulate(df_variant, headers="keys", tablefmt="github"))
        
        # save summary for each variant
        summary_path = f"../results/model_selection/baseline_vs_{variant_name}_summary_{RESULT}.csv"
        df_variant.to_csv(summary_path, index=False)
        print(f"\nSaved summary to: {summary_path}")

    elapsed_mins = (time.time() - t0) / 60.0
    print(f"\n{'='*60}")
    print(f"All model selection runs completed in {elapsed_mins:.2f} minutes.")
    print(f"{'='*60}")

