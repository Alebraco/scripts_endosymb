"""Shared evaluation for the RF endosymbiont classifiers.

    files/model_<label>.joblib                       (rf, le) tuple
    files/model_<label>_holdout_predictions.csv      predictions
    files/model_<label>_holdout_metrics.json         metrics + importances

"""
import json
import os
import sys

import pandas as pd
from joblib import dump
from sklearn.inspection import permutation_importance
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    recall_score,
    roc_auc_score,
)

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from utils import files_dir


def evaluate_and_save(label, rf, le, train_df, holdout_df, features, rf_params):
    """Evaluate a fitted RF on the holdout"""
    y_holdout = le.transform(holdout_df['Group'])
    X_holdout = holdout_df[features].values

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
    oob_error = (1.0 - rf.oob_score_) if hasattr(rf, 'oob_score_') else None
    cm = confusion_matrix(y_holdout, pred)

    gini = dict(zip(features, rf.feature_importances_.tolist()))

    perm = permutation_importance(rf, X_holdout, y_holdout,
                                  n_repeats=30, random_state=28, n_jobs=-1)
    perm_imp = {
        f: {'mean': float(m), 'std': float(s)}
        for f, m, s in zip(features, perm.importances_mean, perm.importances_std)
    }

    true_labels = le.inverse_transform(y_holdout)
    pred_labels = le.inverse_transform(pred)
    pred_df = pd.DataFrame({
        'File': holdout_df['File'].values,
        'Species': holdout_df['Species'].values,
        'true_label': true_labels,
        'pred_label': pred_labels,
        'pred_proba_endosymb': proba_endo,
    })

    mis = pred_df[pred_df['true_label'] != pred_df['pred_label']]
    per_genus = (mis.groupby('Species').size()
                 .sort_values(ascending=False).to_dict())

    results = {
        'model': label,
        'features': list(features),
        'rf_params': rf_params,
        'n_train': int(len(train_df)),
        'n_holdout': int(len(holdout_df)),
        'classes': list(le.classes_),
        'positive_class': 'endosymb_only',
        'metrics': {
            'sensitivity': float(sensitivity),
            'specificity': float(specificity),
            'f1': float(f1),
            'roc_auc': float(roc_auc),
            'oob_error': float(oob_error) if oob_error is not None else None,
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

    results_path = os.path.join(files_dir, 'model_results')
    os.makedirs(results_path, exist_ok=True)
    model_path = os.path.join(results_path, f'model_{label}.joblib')
    pred_path = os.path.join(results_path, f'model_{label}_holdout_predictions.csv')
    metrics_path = os.path.join(results_path, f'model_{label}_holdout_metrics.json')

    dump((rf, le), model_path)
    pred_df.to_csv(pred_path, index=False)
    with open(metrics_path, 'w') as fh:
        json.dump(results, fh, indent=2)

    print('\n' + json.dumps(results['metrics'], indent=2))
    print(f'\nWrote {model_path}')
    print(f'Wrote {pred_path}')
    print(f'Wrote {metrics_path}')

    return results
