"""Create a stratified 90/10 train/holdout split of the RF training set.

Stratification is per (Group, Species) bucket. Buckets with fewer than 3 rows
stay entirely in the training set. Seed is fixed (42) for reproducibility.

Default input path mirrors prediction/models.py:
    endosymb+relatives/feature_files/combined_features.csv
"""
import argparse
import os
import sys

import pandas as pd
from sklearn.model_selection import train_test_split

SEED = 42
TEST_FRAC = 0.10
MIN_BUCKET = 3

DEFAULT_INPUT = os.path.join(
    'endosymb+relatives', 'feature_files', 'combined_features.csv'
)
DEFAULT_OUTDIR = 'files'


def split(df: pd.DataFrame):
    train_parts, holdout_parts = [], []
    rows = []
    for (group, species), bucket in df.groupby(['Group', 'Species'], sort=False):
        n = len(bucket)
        if n < MIN_BUCKET:
            tr, ho = bucket, bucket.iloc[0:0]
        else:
            tr, ho = train_test_split(
                bucket, test_size=TEST_FRAC, random_state=SEED
            )
        train_parts.append(tr)
        if len(ho):
            holdout_parts.append(ho)
        rows.append((group, species, n, len(tr), len(ho)))

    train = pd.concat(train_parts).sort_index()
    holdout = (
        pd.concat(holdout_parts).sort_index()
        if holdout_parts else df.iloc[0:0].copy()
    )
    summary = pd.DataFrame(
        rows, columns=['Group', 'Species', 'n_total', 'n_train', 'n_holdout']
    )
    return train, holdout, summary


def print_summary(df, train, holdout, summary):
    print(f"Seed: {SEED}  test_frac: {TEST_FRAC}  min_bucket: {MIN_BUCKET}")
    print(f"Input rows: {len(df)}  train: {len(train)}  holdout: {len(holdout)}")
    assert len(train) + len(holdout) == len(df), "row count mismatch"

    print("\nPer-class totals:")
    for g in sorted(df['Group'].unique()):
        n = (df['Group'] == g).sum()
        nt = (train['Group'] == g).sum()
        nh = (holdout['Group'] == g).sum()
        frac = nh / n if n else 0.0
        print(f"  {g}: total={n}  train={nt}  holdout={nh}  holdout_frac={frac:.3f}")

    print("\nPer-(Group, Species) split:")
    with pd.option_context('display.max_rows', None, 'display.width', 120):
        print(summary.to_string(index=False))

    held_out_genera = (summary['n_holdout'] > 0).sum()
    held_in_train_only = ((summary['n_holdout'] == 0) & (summary['n_total'] < MIN_BUCKET)).sum()
    print(
        f"\nGenera buckets: {len(summary)}  "
        f"split: {held_out_genera}  training-only (<{MIN_BUCKET}): {held_in_train_only}"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', default=DEFAULT_INPUT, help='Path to combined_features.csv')
    parser.add_argument('--outdir', default=DEFAULT_OUTDIR, help='Output directory')
    args = parser.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"Input CSV not found: {args.input}")

    df = pd.read_csv(args.input)
    required = {'Group', 'Species', 'File'}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"Input CSV missing required columns: {missing}")

    train, holdout, summary = split(df)

    os.makedirs(args.outdir, exist_ok=True)
    train_path = os.path.join(args.outdir, 'training_set.csv')
    holdout_path = os.path.join(args.outdir, 'holdout_locked.csv')
    train.to_csv(train_path, index=False)
    holdout.to_csv(holdout_path, index=False)

    print_summary(df, train, holdout, summary)
    print(f"\nWrote {train_path}")
    print(f"Wrote {holdout_path}  (LOCKED)")


if __name__ == '__main__':
    main()
