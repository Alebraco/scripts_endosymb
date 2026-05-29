"""Train and evaluate Model B (5-feature RF) on the locked train/holdout split.

Inputs:
    files/training_set_B.csv
    files/holdout_locked.csv
Outputs:
    models/model_B.joblib
    results/model_B_holdout_metrics.json
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from joblib import dump
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    recall_score,
    roc_auc_score,
)
from sklearn.preprocessing import LabelEncoder

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from utils import files_dir

MODEL_B_FEATURES = ['Delta_GC2_4', 'GC4', 'AV_Bias', 'Rest_Bias', 'Mean_IGS_Size']

RF_PARAMS = dict(
    n_estimators=500,
    max_depth=None,
    class_weight='balanced',
    random_state=28,
    n_jobs=-1,
    oob_score=True,
)

ROOT = os.path.join(os.path.dirname(__file__), '..')
TRAIN_PATH = os.path.join(files_dir, 'training_set_B.csv')
HOLDOUT_PATH = os.path.join(files_dir, 'holdout_locked.csv')
MODEL_PATH = os.path.join(ROOT, 'models', 'model_B.joblib')
RESULTS_PATH = os.path.join(ROOT, 'results', 'model_B_holdout_metrics.json')


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
    y_holdout = le.transform(holdout['Group'])

    X_train = train[MODEL_B_FEATURES].values
    X_holdout = holdout[MODEL_B_FEATURES].values

    rf = RandomForestClassifier(**RF_PARAMS)
    rf.fit(X_train, y_train)

    endo_idx = list(le.classes_).index('endosymb_only')
    rel_idx = list(le.classes_).index('relatives_only')

    pred = rf.predict(X_holdout)
    proba_endo = rf.predict_proba(X_holdout)[:, endo_idx]

    sensitivity = recall_score(y_holdout, pred, pos_label=endo_idx, average='binary')
    specificity = recall_score(y_holdout, pred, pos_label=rel_idx, average='binary')
    f1 = f1_score(y_holdout, pred, pos_label=endo_idx, average='binary')
    try:
        roc_auc = roc_auc_score((y_holdout == endo_idx).astype(int), proba_endo)
    except ValueError:
        roc_auc = float('nan')
    accuracy = accuracy_score(y_holdout, pred)
    oob_error = 1.0 - rf.oob_score_
    cm = confusion_matrix(y_holdout, pred)

    gini = dict(zip(MODEL_B_FEATURES, rf.feature_importances_.tolist()))

    perm = permutation_importance(rf, X_holdout, y_holdout,
                                  n_repeats=30, random_state=28, n_jobs=-1)
    perm_imp = {
        f: {'mean': float(m), 'std': float(s)}
        for f, m, s in zip(MODEL_B_FEATURES, perm.importances_mean, perm.importances_std)
    }

    pred_labels = le.inverse_transform(pred)
    true_labels = le.inverse_transform(y_holdout)
    mis = pd.DataFrame({'Species': holdout['Species'].values,
                        'true': true_labels, 'pred': pred_labels})
    mis = mis[mis['true'] != mis['pred']]
    per_genus = (mis.groupby('Species').size()
                 .sort_values(ascending=False).to_dict())

    results = {
        'model': 'B',
        'features': MODEL_B_FEATURES,
        'rf_params': RF_PARAMS,
        'n_train': int(len(train)),
        'n_holdout': int(len(holdout)),
        'classes': list(le.classes_),
        'positive_class': 'endosymb_only',
        'metrics': {
            'sensitivity': float(sensitivity),
            'specificity': float(specificity),
            'f1': float(f1),
            'roc_auc': float(roc_auc),
            'oob_error': float(oob_error),
            'accuracy': float(accuracy),
        },
        'confusion_matrix': {
            'label_order': list(le.classes_),
            'rows_true_cols_pred': cm.tolist(),
        },
        'gini_importance': gini,
        'permutation_importance': perm_imp,
        'per_genus_misclassified': {k: int(v) for k, v in per_genus.items()},
    }

    os.makedirs(os.path.dirname(MODEL_PATH), exist_ok=True)
    os.makedirs(os.path.dirname(RESULTS_PATH), exist_ok=True)
    dump((rf, le), MODEL_PATH)
    with open(RESULTS_PATH, 'w') as fh:
        json.dump(results, fh, indent=2)

    print('\n' + json.dumps(results['metrics'], indent=2))
    print(f'\nWrote {MODEL_PATH}')
    print(f'Wrote {RESULTS_PATH}')


if __name__ == '__main__':
    main()
