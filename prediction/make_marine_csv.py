#!/usr/bin/env python3
"""Turn the marine feature output into files/marine_free_livers.csv (Group=relatives_only).

The marine genomes were run through the existing extractor with
`collect_features.py --infer`, which leaves Group='Unknown'. These are all
free-living relatives, so we set Group='relatives_only' and verify the schema
matches files/training_set.csv exactly before writing.

Inputs:
    <path>/feature_files/combined_features.csv   (default path: marine_free_livers)
    files/training_set.csv                        (canonical header reference)
Output:
    files/marine_free_livers.csv
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from utils import files_dir

GROUP_LABEL = 'relatives_only'


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--path', default='marine_free_livers',
                        help='Base dir holding feature_files/combined_features.csv')
    args = parser.parse_args()

    marine_csv = os.path.join(args.path, 'feature_files', 'combined_features.csv')
    train_csv = os.path.join(files_dir, 'training_set.csv')
    out_csv = os.path.join(files_dir, 'marine_free_livers.csv')

    for p in (marine_csv, train_csv):
        if not os.path.exists(p):
            sys.exit(f'Required input not found: {p}')

    marine = pd.read_csv(marine_csv)
    train_cols = list(pd.read_csv(train_csv, nrows=0).columns)

    if list(marine.columns) != train_cols:
        sys.exit(
            'Schema mismatch between marine features and training_set.csv.\n'
            f'  training_set.csv: {train_cols}\n'
            f'  marine          : {list(marine.columns)}'
        )

    n_unknown = int((marine['Group'] == 'Unknown').sum())
    marine['Group'] = GROUP_LABEL
    print(f'Relabeled {n_unknown} rows from Unknown -> {GROUP_LABEL} '
          f'({len(marine)} marine rows total).')

    os.makedirs(files_dir, exist_ok=True)
    marine.to_csv(out_csv, index=False)
    print(f'Wrote {out_csv}\n')

    # Full dump for GATE 1 sanity check (do the values look free-living-like?).
    with pd.option_context('display.max_rows', None, 'display.width', None,
                           'display.max_columns', None):
        print(marine.to_string(index=False))


if __name__ == '__main__':
    main()
