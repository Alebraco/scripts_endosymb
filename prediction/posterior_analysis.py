import os
import sys
import pandas as pd
import arviz as az
import matplotlib.pyplot as plt
import numpy as np
from sklearn.mixture import BayesianGaussianMixture

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from utils import POSTER_RCPARAMS

# Load
df = pd.read_csv('endosymb+relatives/feature_files/processed/processed_bayesian_data.csv')
trace = az.from_netcdf('traces/full_mvn_trace.nc')

feature_columns = ['GC4','AV_Bias','Rest_Bias',
                   'Transposase_Per_Gene','Mean_IGS_Size']

# Rebuild indices (same as in run script)
df['clade_idx'] = df['Species'].astype('category').cat.codes
mask = df['is_endosymbiont'] == 1

# Extract theta posteriors
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
print(f"Exported {len(result)} genomes")


# Bayesian GMM with Dirichlet Process
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

print(f"\nComponent weights: {np.round(weights, 3)}")
print(f"Effective clusters: {effective_k}")
print(f"Means of effective clusters: "
f"{np.round(sorted(bgmm.means_[effective_mask].flatten()), 3)}")

# Assign labels
endo_rows['bgmm_cluster'] = bgmm.predict(endo_theta)
endo_rows['bgmm_prob'] = bgmm.predict_proba(endo_theta).max(axis=1)


# Compute mean theta for each cluster
cluster_means = endo_rows.groupby('bgmm_cluster')['theta_mean'].mean().sort_values()

stage_names = ['recent','transitional','reduced','ancient'][:effective_k]
label_map = {cluster: stage_names[i] for i, cluster in enumerate(cluster_means.index)
            if weights[cluster] > 0.01}

endo_rows['stage'] = endo_rows['bgmm_cluster'].map(label_map)

print("\nCluster summary:")
print(endo_rows.groupby('stage').agg(
    count=('theta_mean','size'),
    theta_mean=('theta_mean','mean'),
    theta_sd=('theta_sd','mean')
).sort_values('theta_mean'))

print("\n>Species composition by cluster")
for stage in ['recent','transitional','reduced','ancient']:
    subset = endo_rows[endo_rows['stage'] == stage]
    top_species = subset['Species'].value_counts().head(8)
    print(f"\n{stage.upper()} (θ̄={round(subset['theta_mean'].mean(),2)}, n={len(subset)}):")
    print(top_species)

endo_rows.to_csv('endosymbiont_clusters.csv', index=False)

stages_present = endo_rows['stage'].dropna().unique()
colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(stages_present)))
color_map = dict(zip(sorted(stages_present,
            key=lambda s: endo_rows[endo_rows['stage']==s]['theta_mean'].mean()),
            colors))

PLAIN_COLOR = '#2f6690'

feature_labels = {
    'GC4': 'GC4',
    'AV_Bias': 'AV Bias',
    'Rest_Bias': 'Rest Bias',
    'Transposase_Per_Gene': 'Transposase Per Gene',
    'Mean_IGS_Size': 'Mean IGS Size',
}

stage_order = ['recent','transitional','reduced','ancient']
stage_colors = dict(zip(stage_order, plt.cm.viridis(np.linspace(0.2, 0.9, 4))))

endo_only = result[result['is_endosymbiont'] == 1].merge(
    endo_rows[['File','stage']], on='File', how='left'
)

beta_decay_mean = trace.posterior['beta_decay'].mean(dim=['chain','draw']).values
mu_global_mean = trace.posterior['mu_global'].mean(dim=['chain','draw']).values

for cluster_coloring in [True, False]:
    suffix = 'clusters' if cluster_coloring else 'plain'

    # Figure 1: Theta Continuum
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
            plt.tight_layout()
            plt.savefig('fig_D_theta_continuum.png')
            plt.close()

    # Figure 2: Inverted-U
    fig, ax = plt.subplots(figsize=(8, 6))
    if cluster_coloring:
        for stage, color in color_map.items():
            subset = endo_rows[endo_rows['stage'] == stage]
            ax.scatter(subset['theta_mean'], subset['theta_sd'],
                       alpha=0.5, s=20, color=color, label=stage.capitalize())
        ax.legend(fontsize=11)
    else:
        ax.scatter(endo_rows['theta_mean'], endo_rows['theta_sd'],
                   alpha=0.5, s=20, color=PLAIN_COLOR)
    ax.set_xlabel('Posterior Mean θ', fontsize=12)
    ax.set_ylabel('Posterior SD θ', fontsize=12)
    ax.set_title('Posterior Uncertainty vs Transition Score', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'figure2_inverted_u_{suffix}.png', dpi=300)
    plt.close()

    # Figure 3: Feature Decay Profiles
    fig, axes = plt.subplots(1, 5, figsize=(20, 4.5), sharey=False)
    for i, feature in enumerate(feature_columns):
        ax = axes[i]
        if cluster_coloring:
            for stage in stage_order:
                subset = endo_only[endo_only['stage'] == stage]
                ax.scatter(subset['theta_mean'], subset[feature],
                           alpha=0.5, s=15, color=stage_colors[stage],
                           label=stage.capitalize() if i == 0 else None)
        else:
            ax.scatter(endo_only['theta_mean'], endo_only[feature],
                       alpha=0.5, s=15, color=PLAIN_COLOR,
                       label='Endosymbiont' if i == 0 else None)
        x_line = np.linspace(0, 2, 50)
        y_line = mu_global_mean[i] + beta_decay_mean[i] * x_line
        ax.plot(x_line, y_line, 'k--', linewidth=1.5,
                label='Model prediction' if i == 0 else None)
        ax.set_xlabel('θ (Transition Score)', fontsize=11)
        ax.set_ylabel(f'{feature_labels[feature]} (z-scored)', fontsize=11)
        ax.set_title(f'β={beta_decay_mean[i]:.2f}', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
        if feature == 'Mean_IGS_Size':
            ax.set_ylim(-3, 6)
        if feature == 'Transposase_Per_Gene':
            ax.set_ylim(-2, 4)
    axes[0].legend(fontsize=10, loc='best')
    fig.suptitle('Feature Decay Profiles Along the Transition Continuum',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(f'figure3_decay_profiles_{suffix}.png', dpi=300, bbox_inches='tight')
    plt.close()

# Figure 4: beta_decay credible intervals
samples = trace.posterior['beta_decay'].values.reshape(-1, len(feature_columns))
means   = samples.mean(axis=0)
lo95, hi95 = np.quantile(samples, [0.025, 0.975], axis=0)

order = np.argsort(means)  # most negative at top when reversed
fig, ax = plt.subplots(figsize=(8, 4))

COLOR = '#2f6690'
for rank, i in enumerate(order[::-1]):
    y = rank
    ax.plot([lo95[i], hi95[i]], [y, y], color=COLOR, lw=2.0,
            marker='|', markersize=8, markeredgewidth=2.0)
    ax.scatter(means[i], y, color=COLOR, s=60, zorder=5, edgecolor='white', lw=1.5)

ax.set_axisbelow(True)
ax.xaxis.grid(True, color='#dde3ea', linewidth=0.8, zorder=0)
ax.axvline(0, color='#555f6b', lw=1.5, linestyle='--', zorder=1)
ax.set_yticks(range(len(feature_columns)))
ax.set_yticklabels([feature_labels[feature_columns[i]] for i in order[::-1]], fontweight='bold')
ax.set_xlabel(r'$\beta_\mathrm{decay}$ (posterior, z-score units)')
ax.set_title('Feature Decay Effects', fontweight='bold', fontsize=14)
for s in ('top', 'right', 'left'):
    ax.spines[s].set_visible(False)
plt.tight_layout()
plt.savefig('figure4_beta_decay_ci.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved figure4_beta_decay_ci.png")