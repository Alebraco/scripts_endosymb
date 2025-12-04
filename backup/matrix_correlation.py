#!/usr/bin/env python3
from skbio.stats.distance import mantel

def matrix_correlation(x_df, y_df):
    '''Computes correlation between two individual delta matrices
    args:
        x_df, y_df: delta matrix (one species only)
    returns:
        (correlation, p_value, n) or None if skipped
    '''
    # Retrieve the order of the IDs
    ordered_ids = x_df.index.tolist()
    y_df = y_df.reindex(index=ordered_ids, columns=ordered_ids)

    if len(x_df) < 3:
        return None

    correlation, p_value, n = mantel(x_df, y_df, method='pearson', permutations=10000)
    p_value_rounded = round(p_value, 4)
    return correlation, p_value_rounded, n