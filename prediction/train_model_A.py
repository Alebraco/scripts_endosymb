"""Train and evaluate Model A (6-feature RF) on the locked train/holdout split.

Inputs:
    files/training_set.csv
    files/holdout_locked.csv
Outputs (via prediction/_model_eval.py):
    files/model_A.joblib
    files/model_A_holdout_predictions.csv
    files/model_A_holdout_metrics.json
"""
import os
import sys

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from utils import feature_columns, files_dir
from model_eval import evaluate_and_save

RF_PARAMS = dict(
    n_estimators=500,
    max_depth=None,
    class_weight='balanced',
    random_state=28,
    n_jobs=-1,
    oob_score=True,
)

TRAIN_PATH = os.path.join(files_dir, 'training_set.csv')
HOLDOUT_PATH = os.path.join(files_dir, 'holdout_locked.csv')


def main():
    for p in (TRAIN_PATH, HOLDOUT_PATH):
        if not os.path.exists(p):
            sys.exit(f'Required input not found: {p}')

    train = pd.read_csv(TRAIN_PATH)
    holdout = pd.read_csv(HOLDOUT_PATH)
    print(f'Training rows: {len(train)}  holdout rows: {len(holdout)}')
    print(f'Features: {feature_columns}')

    le = LabelEncoder()
    le.fit(pd.concat([train['Group'], holdout['Group']], ignore_index=True))
    y_train = le.transform(train['Group'])
    X_train = train[feature_columns].values

    rf = RandomForestClassifier(**RF_PARAMS)
    rf.fit(X_train, y_train)

    evaluate_and_save('A', rf, le, train, holdout, feature_columns, RF_PARAMS)


if __name__ == '__main__':
    main()
