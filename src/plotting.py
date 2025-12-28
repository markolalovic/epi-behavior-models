# src/plotting.py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from pyabc.visualization import plot_kde_matrix

from config import PLOT_LABELS, EXAMPLE_LOCATION
from create_fits_data import generate_predictive_trajectories

def setup_plotting(font_size=11):
    plt.rcParams.update({
        "text.usetex": True,
        "font.size": font_size,
        "axes.labelsize": font_size,
        "xtick.labelsize": font_size - 5,
        "ytick.labelsize": font_size - 5
    })

def plot_posteriors(df_post, w_post, truth, label, 
                    rangex=40, rangey=40):
    """Plots weighted KDE matrix of marginals and joint posteriors.
    """
    n_par = df_post.shape[1]
    fig, axes = plt.subplots(n_par, n_par, figsize=(n_par*2, n_par*2))
    
    plot_kde_matrix(
        df_post, w_post, 
        arr_ax=axes, 
        refval=truth, 
        refval_color='red', 
        names=PLOT_LABELS, 
        colorbar=False,
        numx=rangex, 
        numy=rangey  
    )
    # clean-up
    for i in range(n_par):
        for j in range(n_par):
            ax = axes[i, j]
            if j > i:
                ax.set_visible(False) # lower triangular matrix only
            else:
                if i > j:
                    ax.set_rasterized(True) # rasterize the density polygons
                # format ticks
                if i == n_par - 1:
                    # increased precision for zeta ticks, e.g. 0.012
                    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
                if j == 0:
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    out_path = f"../figures/examples/posteriors_{label}.svg"
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\nPosteriors figure saved to: {out_path}")

def plot_fits(df_post, w_post, model_cfg, obs_deaths, label):
    '''Plots "fit figure" by sampling estimated posterior: 
      - weighted median trajectory 
      - and bands around it (pointwise uncertainty intervals) 
    '''
    post_df_full = df_post.copy()
    post_df_full["weight"] = w_post
    
    # generate ensemble of trajectories
    trajs = generate_predictive_trajectories(
        posterior_df=post_df_full, 
        model_dict=model_cfg, 
        location=EXAMPLE_LOCATION,
        obs_length=len(obs_deaths), 
        n_samples=1000
    )

    q_low, q_med, q_hi = np.quantile(trajs, [0.05, 0.5, 0.95], axis=0)
    days = np.arange(len(obs_deaths))

    plt.figure(figsize=(6, 4))
    plt.fill_between(days, q_low, q_hi, color=model_cfg["color"], alpha=0.3, linewidth=0)
    plt.plot(days, q_med, color=model_cfg["color"], linewidth=1.5)
    plt.scatter(days, obs_deaths, color='black', s=5, alpha=0.6, zorder=3)
    
    plt.xlabel("Day")
    plt.ylabel("Daily Deaths")
    plt.gca().spines[['top', 'right']].set_visible(False)
    plt.grid(False)

    out_path = f"../figures/examples/fit_{label}.svg"
    plt.savefig(out_path, bbox_inches='tight')
    plt.close()
    print(f"\nFits figure saved to: {out_path}")
