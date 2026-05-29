"""Train and evaluate Model B (5-feature RF) on the locked train/holdout split.

Model B drops Transposase_Per_Gene from Model A's feature set and trains on
training_set_B.csv (four endosymbiont genera removed, marine free-livers added).
RF hyperparameters are identical to Model A; random_state=28.

Inputs:
    files/training_set_B.csv
    files/holdout_locked.csv
Outputs (via prediction/_model_eval.py):
    files/model_B.joblib
    files/model_B_holdout_predictions.csv
    files/model_B_holdout_metrics.json
"""
import os
import sys

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from utils import files_dir
from _model_eval import evaluate_and_save

# Model B: 5 features (Transposase_Per_Gene dropped). Local constant on purpose;
# utils.feature_columns is the 6-feature list Model A depends on - leave it alone.
MODEL_B_FEATURES = ['Delta_GC2_4', 'GC4', 'AV_Bias', 'Rest_Bias', 'Mean_IGS_Size']

# Identical to Model A's RF_PARAMS.
RF_PARAMS = dict(
    n_estimators=500,
    max_depth=None,
    class_weight='balanced',
    random_state=28,
    n_jobs=-1,
    oob_score=True,
)

TRAIN_PATH = os.path.join(files_dir, 'training_set_B.csv')
HOLDOUT_PATH = os.path.join(files_dir, 'holdout_locked.csv')


def main():
    for p in (TRAIN_PATH, HOLDOUT_PATH):
        if not os.path.exists(p):
            sys.exit(f'Required input not found: {p}')

    train = pd.read_csv(TRAIN_PATH)
    holdout = pd.read_csv(HOLDOUT_PATH)
    print(f'Training rows: {len(train)}  holdout rows: {len(holdout)}')
    print(f'Features: {MODEL_B_FEATURES}')

    le = LabelEncoder()
    le.fit(pd.concat([train['Group'], holdout['Group']], ignore_index=True))
    y_train = le.transform(train['Group'])
    X_train = train[MODEL_B_FEATURES].values

    rf = RandomForestClassifier(**RF_PARAMS)
    rf.fit(X_train, y_train)

    evaluate_and_save('B', rf, le, train, holdout, MODEL_B_FEATURES, RF_PARAMS)


if __name__ == '__main__':
    main()
