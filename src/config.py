#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/config.py
Configuration file specifies:
- Locations, names, population sizes, 
- Model names, properties for analysis and plotting,
- Fixed model parameters,
- Priors for inferred parameters by ABC-SMC.
"""

from pyabc import Distribution, RV
import numpy as np

# NOTE: temporary
# NOTE: temporary — run model selection only for locations missing
# from results/model_selection/baseline_vs_*_summary_1.csv.
# LOCATIONS = [
#     'AK', 'AL', 'AR', 'AZ', 'CA', 'FL',
#     'HI', 'ID', 'ME', 'MT', 'ND', 'OR',
#     'SC', 'SD', 'TN', 'TX', 'UT', 'VT',
#     'WV', 'WY'
# ]

# selected US locations
# LOCATIONS = [
#     'CO', 'CT', 'DC', 'DE', 'IA', 
#     'IL', 'IN', 'LA', 'MA', 'MD', 
#     'MI', 'MN', 'NJ', 'NM', 'NY', 
#     'OH', 'PA', 'RI', 'VA', 'WA']

# excluded locations
# LOCATIONS_BAD = [
#     'AK', 'AL', 'AR', 'AZ', 'CA', 'FL',
#     'GA', 'HI', 'ID', 'KS', 'KY', 'ME', 
#     'MO', 'MS', 'MT', 'NC', 'ND', 'NE', 
#     'NH', 'NV', 'OK', 'OR', 'SC', 'SD', 
#     'TN', 'TX', 'UT', 'VT', 'WI', 'WV', 
#     'WY']

LOCATIONS = [
    'AK', 'AL', 'AR', 'AZ', 'CA', 'CO', 'CT', 'DC', 'DE', 'FL',
    'GA', 'HI', 'IA', 'ID', 'IL', 'IN', 'KS', 'KY', 'LA', 'MA',
    'MD', 'ME', 'MI', 'MN', 'MO', 'MS', 'MT', 'NC', 'ND', 'NE',
    'NH', 'NJ', 'NM', 'NV', 'NY', 'OH', 'OK', 'OR', 'PA', 'RI',
    'SC', 'SD', 'TN', 'TX', 'UT', 'VA', 'VT', 'WA', 'WI', 'WV',
    'WY'
]

LOCATIONS_BAD = []

# mapping abbreviations to full names 
LOCATION_NAME = {
    'AK': 'Alaska',
    'AL': 'Alabama',
    'AR': 'Arkansas',
    'AZ': 'Arizona',
    'CA': 'California',
    'CO': 'Colorado',
    'CT': 'Connecticut',
    'DC': 'District of Columbia',
    'DE': 'Delaware',
    'FL': 'Florida',
    'GA': 'Georgia',
    'HI': 'Hawaii',
    'IA': 'Iowa',
    'ID': 'Idaho',
    'IL': 'Illinois',
    'IN': 'Indiana',
    'KS': 'Kansas',
    'KY': 'Kentucky',
    'LA': 'Louisiana',
    'MA': 'Massachusetts',
    'MD': 'Maryland',
    'ME': 'Maine',
    'MI': 'Michigan',
    'MN': 'Minnesota',
    'MO': 'Missouri',
    'MS': 'Mississippi',
    'MT': 'Montana',
    'NC': 'North Carolina',
    'ND': 'North Dakota',
    'NE': 'Nebraska',
    'NH': 'New Hampshire',
    'NJ': 'New Jersey',
    'NM': 'New Mexico',
    'NV': 'Nevada',
    'NY': 'New York',
    'OH': 'Ohio',
    'OK': 'Oklahoma',
    'OR': 'Oregon',
    'PA': 'Pennsylvania',
    'RI': 'Rhode Island',
    'SC': 'South Carolina',
    'SD': 'South Dakota',
    'TN': 'Tennessee',
    'TX': 'Texas',
    'UT': 'Utah',
    'VA': 'Virginia',
    'VT': 'Vermont',
    'WA': 'Washington',
    'WI': 'Wisconsin',
    'WV': 'West Virginia',
    'WY': 'Wyoming',
}

# population size from JHU CSSE 
POPULATION_SIZE = {
    'AK': 740995,
    'AL': 4903185,
    'AR': 3017804,
    'AZ': 7278717,
    'CA': 39512223,
    'CO': 5758736,
    'CT': 3565287,
    'DC': 705749,
    'DE': 973764,
    'FL': 21477737,
    'GA': 10617423,
    'HI': 1415872,
    'IA': 3155070,
    'ID': 1787065,
    'IL': 12671821,
    'IN': 6732219,
    'KS': 2913314,
    'KY': 4467673,
    'LA': 4648794,
    'MA': 6892503,
    'MD': 6045680,
    'ME': 1344212,
    'MI': 9986857,
    'MN': 5639632,
    'MO': 6626371,
    'MS': 2976149,
    'MT': 1068778,
    'NC': 10488084,
    'ND': 762062,
    'NE': 1934408,
    'NH': 1359711,
    'NJ': 8882190,
    'NM': 2096829,
    'NV': 3080156,
    'NY': 19453561,
    'OH': 11689100,
    'OK': 3956971,
    'OR': 4217737,
    'PA': 12801989,
    'RI': 1059361,
    'SC': 5148714,
    'SD': 884659,
    'TN': 6829174,
    'TX': 28995881,
    'UT': 3205958,
    'VA': 8535519,
    'VT': 623989,
    'WA': 7614893,
    'WI': 5822434,
    'WV': 1792147,
    'WY': 578759,
}

START_DATE = '2020-03-01'
END_DATE   = '2020-07-01'

MODELS = [
    {
        'name': 'baseline',
        'display_name': 'Baseline',
        'beta_form': 'constant',
        'is_behavioral': False,
        'color': '#0072B2',
    },
    {
        'name': 'behavior_mixed',
        'display_name': 'Behavioral (Mixed)',
        'beta_form': 'mixed',
        'is_behavioral': True,
        'color': '#E69F00',
    },
    {
        'name': 'behavior_exp',
        'display_name': 'Behavioral (Exponential)',
        'beta_form': 'exp',
        'is_behavioral': True,
        'color': '#D55E00',
    },
    {
        'name': 'behavior_rational',
        'display_name': 'Behavioral (Rational)',
        'beta_form': 'rational',
        'is_behavioral': True,
        'color': '#CC79A7',
    },
]

FIXED_PARAMS = {
    'N': None,     # use POPULATION_SIZE[LOCATION]
    'E0': 1.0,
    'R0': 0.0,
    'D0': 0.0,
    'sigma': 1/3,  # mean latent period of 3 days
    'gamma': 1/10, # mean infectious period of 10 days
}

PRIORS = Distribution(
    theta_pi0 = RV('uniform', np.log(1e-8), np.log(1e-3) - np.log(1e-8)),
    R0        = RV('uniform', 1.2, 6.0 - 1.2),
    delta     = RV('uniform', np.log(1e-6), np.log(1e-2) - np.log(1e-6)),
    zeta      = RV('uniform', 0.0, 0.05),
)

# ABC-SMC hyperparameters
ABCSMC_CONFIG = {
    'population_size': 1000,
    'quantile_alpha': 0.3,
    'transition': 'MultivariateNormalTransition',
    'sampler': 'MulticoreEvalParallelSampler',
}
# # NOTE: for testing
# ABCSMC_CONFIG = {
#     'population_size': 50,
#     'quantile_alpha': 0.3,
#     'transition': 'MultivariateNormalTransition',
#     'sampler': 'MulticoreEvalParallelSampler',
# }


# stopping criteria
ABC_RUN_CONFIG = {
    'max_nr_populations': 10,
    'max_walltime_s': 30 * 60,
    'max_total_nr_simulations': 120000,
    'min_eps_diff': 1e-3
}
# # NOTE: for testing
# ABC_RUN_CONFIG = {
#     'max_nr_populations': 2,
#     'max_walltime_s': 5 * 60,
#     'max_total_nr_simulations': 2000,
#     'min_eps_diff': 1e-3
# }

# ABC-SMC model selection hyperparameters 
MODEL_SELECTION_CONFIG = {
    'population_size': 3000,
    'max_nr_populations': 15,
    'quantile_alpha': 0.5,
    'stab_delta': 0.02, # change threshold for g^*
    'stab_k': 3         # consecutive generations for g^*
}

# synthetic experiment settings
SYNTHETIC_CONFIG = {
    'ground_truth_location': 'MA',
    'zeta_grid': [0.001, 0.005, 0.009, 0.012, 0.016, 0.02],
    'noise_sd': 0.025
}

# example settings
# plausible ground-truth values
# see extract_parameters from utils.py
EXAMPLE_LOCATION = 'IN'
EXAMPLE_N = 6732219 # IN, Indiana
EXAMPLE_PARAMS = {
    'theta_pi0': -8.208551,
    'R0': 4.069813,
    'delta': -9.913254,
    'zeta': 0.012352,
}

# plot rendering settings
PLOT_LABELS = {
    "R0": r"$R_0$",
    "theta_pi0": r"$\log(\pi_0)$",
    "delta": r"$\log(\delta)$",
    "zeta": r"$\zeta$"
}

# RESULTS encoding
# NOTE: temporary
# RESULT = 2

RESULT = 1