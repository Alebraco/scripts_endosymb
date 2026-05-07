"""Train Model A (original 6-feature RF) on the locked train/holdout split.

Inputs:
    files/training_set.csv   (from make_holdout_split.py)
    files/holdout_locked.csv (from make_holdout_split.py)

Both inputs are already filtered to Num_Contigs <= 50 upstream by
collect_features.py.

Outputs:
    files/model_A.pkl
    files/model_A_holdout_predictions.csv
    files/model_A_metrics.txt
"""
import os
import sys

import numpy as np
import pandas as pd
from joblib import dump
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.preprocessing import LabelEncoder

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from utils import feature_columns, files_dir

RF_PARAMS = dict(
    n_estimators=500,
    max_depth=None,
    class_weight='balanced',
    random_state=28,
    n_jobs=-1,
)

TRAIN_PATH = os.path.join(files_dir, 'training_set.csv')
HOLDOUT_PATH = os.path.join(files_dir, 'holdout_locked.csv')
MODEL_PATH = os.path.join(files_dir, 'model_A.pkl')
PRED_PATH = os.path.join(files_dir, 'model_A_holdout_predictions.csv')
METRICS_PATH = os.path.join(files_dir, 'model_A_metrics.txt')


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
    y_holdout = le.transform(holdout['Group'])

    X_train = train[feature_columns].values
    X_holdout = holdout[feature_columns].values

    rf = RandomForestClassifier(**RF_PARAMS)
    rf.fit(X_train, y_train)

    pred = rf.predict(X_holdout)
    proba = rf.predict_proba(X_holdout)

    endo_idx = list(le.classes_).index('endosymb_only')
    proba_endo = proba[:, endo_idx]

    acc = accuracy_score(y_holdout, pred)
    prec = precision_score(y_holdout, pred, average='binary', pos_label=endo_idx)
    rec = recall_score(y_holdout, pred, average='binary', pos_label=endo_idx)
    f1 = f1_score(y_holdout, pred, average='binary', pos_label=endo_idx)
    try:
        auc = roc_auc_score((y_holdout == endo_idx).astype(int), proba_endo)
    except ValueError:
        auc = float('nan')
    cm = confusion_matrix(y_holdout, pred)
    report = classification_report(y_holdout, pred, target_names=le.classes_)

    pred_df = pd.DataFrame({
        'File': holdout['File'].values,
        'Species': holdout['Species'].values,
        'true_label': le.inverse_transform(y_holdout),
        'pred_label': le.inverse_transform(pred),
        'pred_proba_endosymb': proba_endo,
    })
    pred_df.to_csv(PRED_PATH, index=False)

    miscls = pred_df[pred_df['true_label'] != pred_df['pred_label']]
    per_genus = (
        miscls.groupby('Species')
        .size()
        .reset_index(name='n_misclassified')
        .sort_values('n_misclassified', ascending=False)
    )

    lines = []
    lines.append('Model A: 6-feature RF on holdout')
    lines.append('')
    lines.append(f'Hyperparameters: {RF_PARAMS}')
    lines.append(f'Features: {feature_columns}')
    lines.append(f'Train rows: {len(train)}   Holdout rows: {len(holdout)}')
    lines.append(f'Classes: {list(le.classes_)}  '
                 f'(positive class for binary metrics: endosymb_only)')
    lines.append('')
    lines.append(f'Accuracy : {acc:.4f}')
    lines.append(f'Precision: {prec:.4f}')
    lines.append(f'Recall   : {rec:.4f}')
    lines.append(f'F1       : {f1:.4f}')
    lines.append(f'ROC-AUC  : {auc:.4f}')
    lines.append('')
    lines.append('Confusion matrix (rows=true, cols=pred, label order = classes above):')
    lines.append(np.array2string(cm))
    lines.append('')
    lines.append('Classification report:')
    lines.append(report)
    lines.append('')
    lines.append('Per-genus misclassifications:')
    if per_genus.empty:
        lines.append('  none')
    else:
        lines.append(per_genus.to_string(index=False))

    metrics_text = '\n'.join(lines)
    with open(METRICS_PATH, 'w') as fh:
        fh.write(metrics_text + '\n')
    print(metrics_text)

    dump((rf, le), MODEL_PATH)
    print(f'\nWrote {MODEL_PATH}')
    print(f'Wrote {PRED_PATH}')
    print(f'Wrote {METRICS_PATH}')


if __name__ == '__main__':
    main()
