"""models.py
Model definitions, run model function, and compute_final_size_until_extinction.
"""

import numpy as np
from scipy.integrate import odeint
from utils import trailing_ma7

def seird_model(y, t, beta0, delta, sigma, gamma, zeta, beta_form):
    r"""
    SEIRD ODE with instantaneous behavioral feedback M(t) = delta * I(t).
    beta_form \in {"constant", "exp", "rational", "mixed"}
    """
    S, E, I, R, D = y
    S = max(S, 0.0); E = max(E, 0.0); I = max(I, 0.0); R = max(R, 0.0)

    N = S + E + I + R
    if N <= 0.0:
        return np.array([0.0, 0.0, 0.0, 0.0, 0.0])

    if beta_form == "constant":
        adjusted_beta = beta0
    else:
        M_t = delta * I
        arg = zeta * M_t
        # clip arg to avoid overflow in exp
        arg_exp = np.clip(arg, -50.0, 50.0)
        arg_den = np.clip(arg, -0.99, 50.0)
        if beta_form == "exp":
            adjusted_beta = beta0 * np.exp(-arg_exp)
        elif beta_form == "rational":
            adjusted_beta = beta0 / (1.0 + arg_den)
        elif beta_form == "mixed":
            adjusted_beta = beta0 * np.exp(-arg_exp) / (1 + arg_den)
        else:
            adjusted_beta = beta0

    force = adjusted_beta * S * I / N
    dS = -force
    dE =  force - sigma * E
    dI =  sigma * E - (gamma + delta) * I
    dR =  gamma * I
    dD =  delta * I

    return np.array([dS, dE, dI, dR, dD])

def run_seird_model(params, 
                    fixed_params, 
                    obs_length, 
                    model_dict,
                    smoothing_window_k=7):
    """
    Runner for the SEIRD model.

    Simulates extra `pad = k-1` days to absorb trailing-MA(k) padding,
    smooths with trailing_ma7 (k=7), then drops the first `pad` values
    so that the returned series has length `obs_length`.

    NOTE: `smoothing_window_k` is consistent with trailing_ma7 window 7
    """
    beta_form     = model_dict['beta_form'] # "constant", "exp", or "rational"
    is_behavioral = model_dict['is_behavioral']

    # fixed parameters
    sigma = float(fixed_params["sigma"])
    gamma = float(fixed_params["gamma"])
    N     = float(fixed_params["N"])
    E0    = float(fixed_params["E0"])
    R0_0  = float(fixed_params["R0"])
    D0    = float(fixed_params["D0"])

    # inferred parameters
    pi0      = float(np.exp(params["theta_pi0"])) # log prevalence -> prevalence
    delta    = float(np.exp(params["delta"]))
    R0_param = float(params["R0"])
    zeta     = float(params["zeta"]) if is_behavioral else 0.0

    # derived
    beta0 = R0_param * (gamma + delta)
    I0    = max(1.0, N * pi0)
    S0    = N - (E0 + I0 + R0_0)
    if S0 < 0.0:
        raise ValueError(f"S0 is negative: {S0}.")
    
    initial_state = [S0, E0, I0, R0_0, D0]


    # --- simulate extra days to absorb trailing-MA padding ---
    # pad = k - 1 for trailing MA(k)
    pad = smoothing_window_k - 1
    T_sim  = int(obs_length + pad)

    # integrate on dayly grid 0, T_sim inclusive (T_sim+1 points -> T_sim daily diffs)
    model_times = np.arange(0, T_sim + 1, dtype=float)
    ode_args    = (beta0, delta, sigma, gamma, zeta, beta_form)

    result, infodict = odeint(
        seird_model, initial_state, model_times, args=ode_args, full_output=True
    )

    if infodict.get('message', '') != 'Integration successful.' or not np.all(np.isfinite(result)):
        bad = np.full(obs_length, np.inf)
        return {"data": bad, "full_trajectory": result}

    # daily deaths from cumulative D (column 4)
    sim_daily_deaths_raw = np.diff(result[:, 4])  # length T_sim

    # smooth by trailing MA(7): trailing_ma7 with fixed-k = 7
    sim_smoothed_ext = trailing_ma7(sim_daily_deaths_raw)

    # guard against tiny negatives
    sim_smoothed_ext = np.maximum(sim_smoothed_ext, 0.0)

    # drop the leading pad values to align to analysis window
    sim_out = sim_smoothed_ext[pad:]

    # final trim to obs_length
    sim_out = np.asarray(sim_out[:obs_length], dtype=float)

    return {"data": sim_out, "full_trajectory": result}

def compute_final_size_until_extinction(params, fixed_params, model_dict,
                                        total_days=720, i_threshold=1.0):
    """
    Integrate the SEIRD model until I(t) < i_threshold
    (or total_days is reached). Returns final size and t_end.
    """
    sim = run_seird_model(
        params=params,
        fixed_params=fixed_params,
        obs_length=total_days,
        model_dict=model_dict
    )

    traj = sim["full_trajectory"]
    S, I = traj[:, 0], traj[:, 2]
    N = float(fixed_params["N"])

    idx = np.where(I < i_threshold)[0]
    t_end = int(idx[0]) if len(idx) > 0 else len(I) - 1

    S_t_end = float(S[t_end])
    final_size = 1.0 - S_t_end / N

    return t_end, S_t_end, final_size
