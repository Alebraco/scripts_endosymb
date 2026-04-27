#!/usr/bin/env python3
"""
Correlate theta estimate (posterior mean) with:
  1. Median patristic distance to free-living relatives.
  2. Genome size.

Endosymbionts only, since free-living genomes have theta ~ 0 by design.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

from distance_matrix import distance_matrix
from gcsize_dict import genome_gcsize
from utils import (
    load_or_compute,
    load_or_compute_pickle,
    files_dir,
    genome_gcsize_json_path,
)

THETA_CSV = os.path.join("files", "theta_all_genomes.csv")
GROUP = "endosymb+relatives"
OUT_CSV = os.path.join(files_dir, "theta_correlations.csv")
OUT_STATS = os.path.join(files_dir, "theta_correlations_stats.txt")
PLOT_DIR = os.path.join("plots", "theta_correlations")


def resolve_endo_id(file_value, matrix_index):
    """Return the matrix-index entry that matches the theta-table File value.

    The matrix index uses raw accession IDs (relatives carry a '_genomic'
    suffix, endosymbionts do not).
    """

    if file_value in matrix_index:
        return file_value
    candidates = [
        idx for idx in matrix_index
        if "_genomic" not in idx and (idx in file_value or file_value in idx)
    ]
    if len(candidates) == 1:
        return candidates[0]
    return None


def build_table():
    if not os.path.exists(THETA_CSV):
        sys.exit(
            f"Missing {THETA_CSV}. Run prediction/posterior_analysis.py first."
        )

    theta_df = pd.read_csv(THETA_CSV)
    endo = theta_df[theta_df["is_endosymbiont"] == 1].copy()
    print(f"Loaded {len(endo)} endosymbiont rows from {THETA_CSV}")

    distance_pkl = os.path.join(files_dir, f"distances_{GROUP}_protein.pkl")
    distance_matrices = load_or_compute_pickle(
        distance_pkl, distance_matrix, GROUP, mode="protein"
    )
    genome_dataset = load_or_compute(
        genome_gcsize_json_path(GROUP), genome_gcsize, GROUP
    )

    rows = []
    skipped = {"no_species_matrix": 0, "no_relatives": 0,
               "id_unresolved": 0, "no_size_entry": 0}

    for _, row in endo.iterrows():
        species = row["Species"]
        file_value = row["File"]

        matrix = distance_matrices.get(species)
        if matrix is None:
            skipped["no_species_matrix"] += 1
            continue

        relatives = [idx for idx in matrix.index if "_genomic" in idx]
        if not relatives:
            skipped["no_relatives"] += 1
            continue

        endo_id = resolve_endo_id(file_value, matrix.index)
        if endo_id is None:
            skipped["id_unresolved"] += 1
            print(f"  unresolved: {species} / {file_value}")
            continue

        median_distance = float(np.median(
            [matrix.loc[endo_id, rel_id] for rel_id in relatives]
        ))

        sp_entry = genome_dataset.get(species, {})
        size_entry = sp_entry.get(endo_id) or sp_entry.get(file_value)
        if size_entry is None:
            skipped["no_size_entry"] += 1
            print(f"  no size: {species} / {endo_id}")
            continue

        rows.append({
            "Species": species,
            "File": file_value,
            "matrix_id": endo_id,
            "theta_mean": row["theta_mean"],
            "theta_sd": row["theta_sd"],
            "median_distance": median_distance,
            "genome_size": int(size_entry["size"]),
        })

    print(f"Retained {len(rows)} / {len(endo)} endosymbionts. "
          f"Skipped: {skipped}")

    df = pd.DataFrame(rows)
    os.makedirs(files_dir, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV}")
    return df


def run_correlations(df):
    pairs = [
        ("theta_mean", "median_distance"),
        ("theta_mean", "genome_size"),
    ]
    lines = [f"n = {len(df)}"]
    results = {}
    for x, y in pairs:
        sub = df[[x, y]].dropna()
        rho, p_s = spearmanr(sub[x], sub[y])
        r, p_p = pearsonr(sub[x], sub[y])
        results[(x, y)] = {"spearman": (rho, p_s), "pearson": (r, p_p),
                           "n": len(sub)}
        lines.append(
            f"{x} vs {y}  (n={len(sub)})  "
            f"Spearman rho={rho:.4f} p={p_s:.3e}  "
            f"Pearson r={r:.4f} p={p_p:.3e}"
        )
    text = "\n".join(lines)
    print(text)
    with open(OUT_STATS, "w") as f:
        f.write(text + "\n")
    print(f"Wrote {OUT_STATS}")
    return results


def make_plots(df, results):
    os.makedirs(PLOT_DIR, exist_ok=True)

    rho, p = results[("theta_mean", "median_distance")]["spearman"]
    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(
        df["theta_mean"], df["median_distance"],
        c=df["genome_size"] / 1e6, cmap="viridis",
        s=40, edgecolors="black", linewidths=0.3, alpha=0.85,
    )
    cbar = plt.colorbar(sc)
    cbar.set_label("Genome size (Mb)", fontsize=12)
    ax.set_xlabel(r"Posterior mean $\theta$", fontsize=13)
    ax.set_ylabel("Median distance to relatives (protein)", fontsize=13)
    ax.set_title(
        rf"$\theta$ vs distance-to-relatives  "
        rf"(Spearman $\rho={rho:.3f}$, p={p:.2e}, n={len(df)})",
        fontsize=13,
    )
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "theta_vs_distance.pdf"))
    plt.close()

    rho, p = results[("theta_mean", "genome_size")]["spearman"]
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(
        df["theta_mean"], df["genome_size"],
        s=40, edgecolors="black", linewidths=0.3, alpha=0.85,
        color="#2f6690",
    )
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _p: f"{x/1e6:.1f} Mb"))
    ax.set_xlabel(r"Posterior mean $\theta$", fontsize=13)
    ax.set_ylabel("Genome size", fontsize=13)
    ax.set_title(
        rf"$\theta$ vs genome size  "
        rf"(Spearman $\rho={rho:.3f}$, p={p:.2e}, n={len(df)})",
        fontsize=13,
    )
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "theta_vs_size.pdf"))
    plt.close()
    print(f"Plots written to {PLOT_DIR}/")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true",
                    help="Recompute the merged table even if cached.")
    args = ap.parse_args()

    if os.path.exists(OUT_CSV) and not args.force:
        print(f"Loading cached {OUT_CSV} (use --force to recompute).")
        df = pd.read_csv(OUT_CSV)
    else:
        df = build_table()

    if df.empty:
        sys.exit("No rows retained, nothing to correlate.")

    results = run_correlations(df)
    make_plots(df, results)


if __name__ == "__main__":
    main()
