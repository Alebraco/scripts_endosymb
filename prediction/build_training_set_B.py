#!/usr/bin/env python3
"""Build files/training_set_B.csv for Model B.

Starting from training_set.csv: drop four endosymbiont genera, then append the
marine free-living relatives. Holdout (holdout_locked.csv) is NOT touched.

Inputs:
    files/training_set.csv
    files/marine_free_livers.csv   (from make_marine_csv.py)
Output:
    files/training_set_B.csv
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from utils import files_dir

# Species to drop, matched case-insensitively as substrings of the Species column.
# 'Serratia symbiotica' is the full phrase so free-living Serratia are kept.
DROP_PATTERNS = ['Sodalis', 'Serratia symbiotica', 'Hamiltonella', 'Arsenophonus']

TRAIN_PATH = os.path.join(files_dir, 'training_set.csv')
MARINE_PATH = os.path.join(files_dir, 'marine_free_livers.csv')
OUT_PATH = os.path.join(files_dir, 'training_set_B.csv')


def main():
    for p in (TRAIN_PATH, MARINE_PATH):
        if not os.path.exists(p):
            sys.exit(f'Required input not found: {p}')

    train = pd.read_csv(TRAIN_PATH)
    print(f'training_set.csv: {len(train)} rows')

    species = train['Species'].fillna('')
    drop_mask = pd.Series(False, index=train.index)
    for pat in DROP_PATTERNS:
        m = species.str.contains(pat, case=False, regex=False)
        print(f'  matched "{pat}": {int(m.sum())} rows')
        drop_mask |= m

    kept = train[~drop_mask].copy()
    print(f'Dropped {int(drop_mask.sum())} rows; {len(kept)} remain after genus removal.')

    marine = pd.read_csv(MARINE_PATH)
    if list(marine.columns) != list(train.columns):
        sys.exit('Schema mismatch between marine_free_livers.csv and training_set.csv.')
    print(f'marine_free_livers.csv: {len(marine)} rows to append')

    combined = pd.concat([kept, marine], ignore_index=True)
    print(f'Final training_set_B.csv: {len(combined)} rows')

    combined.to_csv(OUT_PATH, index=False)
    print(f'Wrote {OUT_PATH}')


if __name__ == '__main__':
    main()
