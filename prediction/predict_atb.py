#!/usr/bin/env python3
"""
Apply the trained Random Forest endosymb/free-living classifier to features
collected from the ATB bacteria dataset.

Inputs:
    --features-csv  Path to combined_features.csv produced by collect_features.py
    --model-path    Path to serialized (rf, scaler, le) tuple from models.py
    --prob-threshold High-confidence cutoff for endosymbiont predictions
                     (default 0.75; matches confidence_threshold in models.py).

Outputs:
    rf_predictions.csv             Per-genome predictions (all rows).
    high_conf_endosymbionts.csv    Filtered list of high-confidence
                                   endosymbiont predictions, with Species.
    high_conf_endosymbiont_species.txt
                                   Unique species names (one per line) that
                                   have at least one high-confidence
                                   endosymbiont call.
"""

import argparse
import os
import sys

import pandas as pd
from joblib import load

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from utils import feature_columns


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--features-csv',
        default='atb_bacteria/feature_files/combined_features.csv',
        help='Path to combined_features.csv',
    )
    parser.add_argument(
        '--model-path',
        default='files/rf_model.joblib',
        help='Path to serialized (rf, scaler, le) joblib',
    )
    parser.add_argument(
        '--prob-threshold',
        type=float,
        default=0.75,
        help='Minimum endosymb probability to flag as high-confidence',
    )
    args = parser.parse_args()

    if not os.path.exists(args.features_csv):
        raise FileNotFoundError(f'Features CSV not found: {args.features_csv}')
    if not os.path.exists(args.model_path):
        raise FileNotFoundError(f'Model file not found: {args.model_path}')

    out_dir = os.path.dirname(os.path.abspath(args.features_csv))

    rf, scaler, le = load(args.model_path)
    df = pd.read_csv(args.features_csv)

    X = df[feature_columns]
    X_scaled = scaler.transform(X)

    probs = rf.predict_proba(X_scaled)
    preds = le.inverse_transform(rf.predict(X_scaled))
    endosymb_idx = list(le.classes_).index('endosymb_only')
    endosymb_prob = probs[:, endosymb_idx]

    results = pd.DataFrame({
        'Species': df['Species'],
        'File': df['File'],
        'Predicted_Label': preds,
        'Endosymb_Probability': endosymb_prob,
        'Confidence': [max(p, 1 - p) for p in endosymb_prob],
        'High_Conf_Endosymb': (
            (preds == 'endosymb_only') & (endosymb_prob >= args.prob_threshold)
        ),
    })

    preds_path = os.path.join(out_dir, 'rf_predictions.csv')
    results.to_csv(preds_path, index=False)
    print(f'Wrote {len(results)} predictions to {preds_path}')

    high_conf = results[results['High_Conf_Endosymb']].copy()
    high_conf = high_conf.sort_values(
        ['Species', 'Endosymb_Probability'], ascending=[True, False]
    )
    high_conf_path = os.path.join(out_dir, 'high_conf_endosymbionts.csv')
    high_conf.to_csv(high_conf_path, index=False)
    print(
        f'Wrote {len(high_conf)} high-confidence endosymbiont predictions '
        f'(p >= {args.prob_threshold}) to {high_conf_path}'
    )

    species_path = os.path.join(out_dir, 'high_conf_endosymbiont_species.txt')
    unique_species = sorted(high_conf['Species'].dropna().unique())
    with open(species_path, 'w') as fh:
        for sp in unique_species:
            fh.write(f'{sp}\n')
    print(
        f'Wrote {len(unique_species)} unique high-confidence endosymbiont '
        f'species to {species_path}'
    )


if __name__ == '__main__':
    main()
