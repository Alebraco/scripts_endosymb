#!/usr/bin/env python3
"""
Project ATB features onto the training PCA space to diagnose prediction quality.

Outputs (in --outdir):
    atb_pca_overlay.pdf        Training clusters + ATB points colored by Endosymb_Probability.
    atb_feature_distributions.pdf  Box plots comparing training vs ATB high-conf predictions.

Usage (on the cluster, from the project root):
    python prediction/visualize_atb_pca.py
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import load
from sklearn.decomposition import PCA

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from utils import feature_columns

BASE = '/rsstu/users/l/ljbobay/recombination/asoneto'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--training-csv',
                   default='endosymb+relatives/feature_files/combined_features.csv')
    p.add_argument('--atb-csv',
                   default=os.path.join(BASE, 'atb_bacteria/feature_files/combined_features.csv'))
    p.add_argument('--predictions-csv',
                   default=os.path.join(BASE, 'atb_bacteria/feature_files/rf_predictions.csv'))
    p.add_argument('--model-path', default='files/rf_model.joblib')
    p.add_argument('--outdir',
                   default=os.path.join(BASE, 'atb_bacteria/feature_files'))
    return p.parse_args()


def load_training(csv_path, scaler):
    df = pd.read_csv(csv_path)
    df['Group'] = df['Group'].replace({
        'endosymb_only': 'Endosymbionts',
        'relatives_only': 'Free-Living Relatives',
    })
    X = pd.DataFrame(scaler.transform(df[feature_columns]), columns=feature_columns)
    return df, X


def plot_pca_overlay(train_df, train_pca, atb_df, atb_pca, var_explained, outdir):
    fig, ax = plt.subplots(figsize=(10, 8))

    palette = {'Endosymbionts': '#FC8D62', 'Free-Living Relatives': '#66C2A5'}
    for group, color in palette.items():
        mask = train_df['Group'].values == group
        ax.scatter(
            train_pca[mask, 0], train_pca[mask, 1],
            c=color, label=f'Training – {group}',
            s=40, alpha=0.35, edgecolors='none',
        )

    probs = atb_df['Endosymb_Probability'].fillna(0.5).values
    sc = ax.scatter(
        atb_pca[:, 0], atb_pca[:, 1],
        c=probs, cmap='RdYlGn_r',
        s=18, alpha=0.7, marker='^',
        vmin=0, vmax=1, label='ATB genomes',
    )
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('Endosymbiont Probability', fontsize=11)

    ax.set_xlabel(f'PC1 ({var_explained[0] * 100:.1f}%)', fontweight='bold')
    ax.set_ylabel(f'PC2 ({var_explained[1] * 100:.1f}%)', fontweight='bold')
    ax.set_title('ATB Predictions Projected onto Training PCA Space', fontweight='bold')
    ax.legend(loc='best', fontsize=9, framealpha=0.7)
    plt.tight_layout()

    out = os.path.join(outdir, 'atb_pca_overlay.pdf')
    fig.savefig(out)
    plt.close(fig)
    print(f'Saved PCA overlay → {out}')


def plot_feature_distributions(train_df, train_X_raw, atb_df, atb_X_raw, outdir):
    high_conf = atb_df['High_Conf_Endosymb'].fillna(False)

    endo_mask = train_df['Group'].values == 'Endosymbionts'
    rel_mask  = train_df['Group'].values == 'Free-Living Relatives'
    hc_mask   = high_conf.values

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes = axes.flatten()

    label_map = {
        'Delta_GC2_4':          'ΔGC2–4',
        'GC4':                  'GC4 Content',
        'AV_Bias':              'Ala/Val Bias',
        'Rest_Bias':            'Rest AA Bias',
        'Transposase_Per_Gene': 'Transposase / Gene',
        'Mean_IGS_Size':        'Mean IGS Size',
    }

    for ax, feat in zip(axes, feature_columns):
        groups = {
            'Training\nEndosymbionts':   train_X_raw[feat].values[endo_mask],
            'Training\nFree-Living':     train_X_raw[feat].values[rel_mask],
            'ATB High-Conf\nPredictions': atb_X_raw[feat].values[hc_mask],
        }
        colors = ['#FC8D62', '#66C2A5', '#8DA0CB']
        positions = range(1, len(groups) + 1)

        bp = ax.boxplot(
            list(groups.values()),
            positions=list(positions),
            patch_artist=True,
            widths=0.5,
            flierprops=dict(marker='.', markersize=3, alpha=0.4),
            medianprops=dict(color='black', linewidth=1.5),
        )
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax.set_xticks(list(positions))
        ax.set_xticklabels(list(groups.keys()), fontsize=8)
        ax.set_title(label_map.get(feat, feat), fontweight='bold')

    plt.suptitle('Feature Distributions: Training vs ATB High-Confidence Predictions',
                 fontweight='bold', y=1.01)
    plt.tight_layout()

    out = os.path.join(outdir, 'atb_feature_distributions.pdf')
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved feature distributions → {out}')


def print_feature_summary(train_df, train_X_raw, atb_df, atb_X_raw):
    high_conf = atb_df['High_Conf_Endosymb'].fillna(False)
    endo_mask = train_df['Group'].values == 'Endosymbionts'
    rel_mask  = train_df['Group'].values == 'Free-Living Relatives'
    hc_mask   = high_conf.values

    def fmt(arr):
        return f'{arr.mean():.3f} ± {arr.std():.3f}  [{arr.min():.3f}, {arr.max():.3f}]'

    col = 32
    print(f'\n{"Feature":<25} {"Train Endosymb (mean±SD)":<{col}} {"Train Free-Living (mean±SD)":<{col}} {"ATB High-Conf (mean±SD)"}')
    print('-' * (25 + col * 2 + 30))
    for feat in feature_columns:
        te = train_X_raw[feat].values[endo_mask]
        fl = train_X_raw[feat].values[rel_mask]
        ah = atb_X_raw[feat].values[hc_mask]
        ah_str = fmt(ah) if hc_mask.any() else 'n/a'
        print(f'{feat:<25} {fmt(te):<{col}} {fmt(fl):<{col}} {ah_str}')

    print(f'\nTotal ATB genomes: {len(atb_df)}')
    print(f'High-confidence endosymbiont predictions: {hc_mask.sum()}')
    if hc_mask.any():
        probs = atb_df.loc[hc_mask, 'Endosymb_Probability']
        print(f'Probability range (high-conf): [{probs.min():.3f}, {probs.max():.3f}], '
              f'median={probs.median():.3f}')


def main():
    args = parse_args()

    for path, label in [(args.training_csv, 'training CSV'),
                         (args.atb_csv, 'ATB features CSV'),
                         (args.predictions_csv, 'predictions CSV'),
                         (args.model_path, 'model file')]:
        if not os.path.exists(path):
            raise FileNotFoundError(f'{label} not found: {path}')

    os.makedirs(args.outdir, exist_ok=True)

    _, scaler, _ = load(args.model_path)

    train_df, X_train_scaled = load_training(args.training_csv, scaler)
    train_X_raw = pd.read_csv(args.training_csv)[feature_columns]
    train_df['Group'] = train_df['Group'].replace({
        'endosymb_only': 'Endosymbionts',
        'relatives_only': 'Free-Living Relatives',
    })

    pca = PCA(n_components=2)
    train_pca = pca.fit_transform(X_train_scaled.values)

    atb_features = pd.read_csv(args.atb_csv)
    preds_df = pd.read_csv(args.predictions_csv)
    atb_df = atb_features.merge(
        preds_df[['File', 'Endosymb_Probability', 'Predicted_Label', 'High_Conf_Endosymb']],
        on='File', how='left',
    )
    atb_X_raw = atb_df[feature_columns]
    atb_X_scaled = scaler.transform(atb_X_raw)
    atb_pca = pca.transform(atb_X_scaled)

    plot_pca_overlay(train_df, train_pca, atb_df, atb_pca,
                     pca.explained_variance_ratio_, args.outdir)
    plot_feature_distributions(train_df, train_X_raw, atb_df, atb_X_raw, args.outdir)
    print_feature_summary(train_df, train_X_raw, atb_df, atb_X_raw)


if __name__ == '__main__':
    main()
