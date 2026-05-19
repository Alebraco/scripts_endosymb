import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.mixture import BayesianGaussianMixture

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from utils import POSTER_RCPARAMS

files_dir = os.path.join(os.path.dirname(__file__), '..', 'files')
clusters_path = os.path.join(files_dir, 'endosymb_clusters.csv')

if os.path.exists(clusters_path):
    endo_rows = pd.read_csv(clusters_path)
else:
    import arviz as az

    df = pd.read_csv('endosymb+relatives/feature_files/processed/processed_bayesian_data.csv')
    trace = az.from_netcdf('traces/full_mvn_trace.nc')

    feature_columns = ['GC4', 'AV_Bias', 'Rest_Bias',
                       'Transposase_Per_Gene', 'Mean_IGS_Size']

    df['clade_idx'] = df['Species'].astype('category').cat.codes
    mask = df['is_endosymbiont'] == 1

    theta_endo_summary = az.summary(trace, var_names=['theta_endo'])
    theta_free_summary = az.summary(trace, var_names=['theta_free'])

    endo_rows = df[mask].copy()
    endo_rows['theta_mean'] = theta_endo_summary['mean'].values
    endo_rows['theta_sd'] = theta_endo_summary['sd'].values

    free_rows = df[~mask].copy()
    free_rows['theta_mean'] = theta_free_summary['mean'].values
    free_rows['theta_sd'] = theta_free_summary['sd'].values

    result = pd.concat([endo_rows, free_rows])
    result.to_csv('theta_all_genomes.csv', index=False)

    endo_theta = endo_rows['theta_mean'].values.reshape(-1, 1)
    bgmm = BayesianGaussianMixture(
        n_components=10,
        weight_concentration_prior_type='dirichlet_process',
        weight_concentration_prior=0.1,
        random_state=42,
        max_iter=500
    )
    bgmm.fit(endo_theta)

    weights = bgmm.weights_
    effective_mask = weights > 0.01
    effective_k = effective_mask.sum()

    endo_rows['bgmm_cluster'] = bgmm.predict(endo_theta)
    endo_rows['bgmm_prob'] = bgmm.predict_proba(endo_theta).max(axis=1)

    cluster_means = endo_rows.groupby('bgmm_cluster')['theta_mean'].mean().sort_values()
    stage_names = ['recent', 'transitional', 'reduced', 'ancient'][:effective_k]
    label_map = {cluster: stage_names[i] for i, cluster in enumerate(cluster_means.index)
                 if weights[cluster] > 0.01}
    endo_rows['stage'] = endo_rows['bgmm_cluster'].map(label_map)

    endo_rows.to_csv(clusters_path, index=False)
    print(f"Saved {clusters_path}")

# Build color map
stages_present = endo_rows['stage'].dropna().unique()
colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(stages_present)))
color_map = dict(zip(
    sorted(stages_present,
           key=lambda s: endo_rows[endo_rows['stage'] == s]['theta_mean'].mean()),
    colors
))

PLAIN_COLOR = '#2f6690'

# Figure 1: Theta Continuum
for cluster_coloring in [True, False]:
    suffix = 'clusters' if cluster_coloring else 'plain'

    fig, ax = plt.subplots(figsize=(10, 5))
    if cluster_coloring:
        for stage, color in color_map.items():
            subset = endo_rows[endo_rows['stage'] == stage]
            ax.hist(subset['theta_mean'], bins=30, alpha=0.7,
                    label=f"{stage.capitalize()} (n={len(subset)})", color=color)
        ax.legend(fontsize=11)
    else:
        ax.hist(endo_rows['theta_mean'], bins=30, alpha=0.7, color=PLAIN_COLOR)
    ax.set_xlabel('Posterior Mean θ (Transition Score)', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Endosymbiont Transition Continuum', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'figure1_theta_continuum_{suffix}.png', dpi=300)
    plt.close()

    if cluster_coloring:
        with plt.rc_context(POSTER_RCPARAMS):
            fig, ax = plt.subplots(figsize=(11, 6))
            for stage, color in color_map.items():
                subset = endo_rows[endo_rows['stage'] == stage]
                ax.hist(subset['theta_mean'], bins=30, alpha=0.7,
                        label=f"{stage.capitalize()} (n={len(subset)})",
                        color=color)
            ax.legend(loc='upper right')
            ax.set_xlabel('Posterior mean θ (transition score)')
            ax.set_ylabel('Count')
            ax.set_title('Endosymbiont Transition Continuum', fontweight='bold')
            plt.tight_layout()
            plt.savefig('fig_D_theta_continuum.png')
            plt.close()
