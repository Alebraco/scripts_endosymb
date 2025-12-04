#!/usr/bin/env python3
import numpy as np
import Bio.Phylo as Phylo
from scipy.cluster.hierarchy import linkage, fcluster
import pandas as pd

# Load tree and calculate distances
tree = Phylo.read("Buchnera.treefile", "newick")
terminals = tree.get_terminals()
dist_matrix = [[tree.distance(t1, t2) for t2 in terminals] for t1 in terminals]
genome_ids = [t.name for t in terminals]

# Cluster genomes (5 clusters by default)
distance_threshold = 0.9
clusters = fcluster(linkage(dist_matrix, method='average'), t=distance_threshold, criterion='distance')

# Save cluster assignments to CSV
results = pd.DataFrame({"Genome_ID": genome_ids, "Cluster": clusters})
results.to_csv("clusters.csv", index=False)
print(f"Saved cluster assignments to clusters.csv")

# Print one representative per cluster
reps = results.groupby("Cluster")["Genome_ID"].first()
print("\nCluster representatives:")
print(reps)

# Calculate average intra-cluster distances
print("\nAverage intra-cluster distances:")
cluster_ids = sorted(results["Cluster"].unique())
for c in cluster_ids:
    idx = results[results["Cluster"] == c].index
    if len(idx) > 1:  # Need at least 2 genomes for distance
        distances = [dist_matrix[i][j] for i in idx for j in idx if i < j]
        avg_dist = np.mean(distances) if distances else 0.0
        print(f"Cluster {c}: {avg_dist:.4f}")
    else:
        print(f"Cluster {c}: 0.0000 (single genome)")
