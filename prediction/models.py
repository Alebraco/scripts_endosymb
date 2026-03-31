import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict
from sklearn.metrics import confusion_matrix, classification_report, accuracy_score, ConfusionMatrixDisplay
from sklearn.inspection import permutation_importance
import plotly.express as px
import numpy as np
from joblib import dump, load
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from utils import feature_columns


feature_dir = 'feature_files'
features_df = 'combined_features.csv'

recent_endosymb = [
    'Arsenophonus',
    'Sodalis',
    'Hamiltonella',
    'Richelia',
    'Serratia symbiotica',
    'Rhizobium',
]

confidence_threshold = 0.75


def is_recent_endosymb(species_name):
    return any(species_name.startswith(r) or r in species_name for r in recent_endosymb)

def plot_correlation(csv_path, outpath):
    plot_path = os.path.join(outpath, 'feature_correlation.pdf')
    df = pd.read_csv(csv_path)

    features = df[feature_columns]

    corr = features.corr()

    plt.figure(figsize=(10, 8))
    sns.heatmap(corr, annot=True, cmap='coolwarm', fmt=".2f", linewidths=0.5)

    plt.title("Feature Correlation Matrix")
    plt.tight_layout()

    plt.savefig(plot_path)
    plt.close()
    
    upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))
    to_drop = [column for column in upper.columns if any(upper[column].abs() > 0.90)]
    if to_drop:
        print(f"Highly correlated features (>0.9): {to_drop}")


def plot_feature_distributions(X, y, outpath):
    plot_df = X.copy()
    plot_df['Group'] = y['Group'].values

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    for i, feature in enumerate(feature_columns):
        ax = axes[i]
        
        sns.boxplot(data=plot_df, x='Group', y=feature, 
                    color='white', width=0.5, fliersize=0, ax=ax)
        
        sns.stripplot(data=plot_df, x='Group', y=feature, 
                      alpha=0.6, jitter=True, hue='Group', legend=False,
                      palette={'Endosymbionts': '#FC8D62', 'Free-Living Relatives': '#66C2A5'}, ax=ax)
        
        ax.set_title(f'{feature}', fontweight='bold')
        ax.set_xlabel('')

    plt.tight_layout()    
    fig.savefig(os.path.join(outpath, f'feature_distributions_grid.pdf'))
    plt.close(fig)
    print("Feature distribution plots saved.")

def run_pca(csv_path, outpath, label_recent=False, interactive=False):
    plot_path = os.path.join(outpath, 'pca_plot.pdf')
    interactive_path = os.path.join(outpath, 'pca_plot.html')
    loadings_path = os.path.join(outpath, 'pca_loadings.csv')
    df = pd.read_csv(csv_path)
    
    df['Group'] = df['Group'].replace({
        'endosymb_only': 'Endosymbionts',
        'relatives_only': 'Free-Living Relatives'
    })
    
    y = ['Group', 'Species', 'File']
    df_ml = df.set_index(y)
    X = df_ml[feature_columns]
    
    # Standardize data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(X_scaled)
    
    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'], index=df_ml.index)
    pca_df = pca_df.reset_index() 
    
    var_explained = pca.explained_variance_ratio_

    if not label_recent:
        plt.figure(figsize=(8, 6))
        sns.scatterplot(
            x='PC1', y='PC2',
            hue='Group',
            palette={'Endosymbionts': '#FC8D62', 'Free-Living Relatives': '#66C2A5'},
            data=pca_df,
            s=50, alpha=0.8, edgecolor='k'
        )
    else:
        suspects = ['Sodalis', 'Serratia symbiotica', 'Richelia', 'Arsenophonus', 'Rickettsia', 'Spiroplasma']
        
        pca_df['Display_Label'] = pca_df.apply(
            lambda row: row['Species'] if row['Species'] in suspects else row['Group'], axis=1
        )

        custom_palette = {
            'Endosymbionts': "#D3D3D3",         
            'Free-Living Relatives': "#D2E1DCC1",  
            'Sodalis': '#E31A1C',                  # Bright Red
            'Serratia symbiotica': "#3C00FF",      # Bright Blue
            'Richelia': '#6A3D9A' ,                 # Deep Purple
            'Arsenophonus': "#B9D70D",                # Bright Orange
            'Rickettsia': "#33A02C",                 # Bright Green
            'Spiroplasma': "#00F7FF",                 # Bright Orange
        }

        bg_mask = pca_df['Display_Label'].isin(['Endosymbionts', 'Free-Living Relatives'])
        plt.figure(figsize=(9, 7))
        sns.scatterplot(
            x='PC1', y='PC2',
            hue='Display_Label',
            palette=custom_palette,
            data=pca_df[bg_mask],
            s=60, alpha=0.5
        )
        sns.scatterplot(
            x='PC1', y='PC2',
            hue='Display_Label',
            palette=custom_palette,
            data=pca_df[~bg_mask],
            style='Group',
            s=60, alpha=0.8
        )

        handles, labels = plt.gca().get_legend_handles_labels()
        clean_labels = []
        clean_handles = []
        for handle, label in zip(handles, labels):
            if label not in ['Group', 'Display_Label'] and label not in clean_labels:
                clean_labels.append(label)
                clean_handles.append(handle)

        plt.legend(clean_handles, clean_labels, loc='best', title='Species/Group')

    plt.title('PCA of Genomic Features', fontsize=16, fontweight='bold')
    plt.xlabel(f'Principal Component 1 ({var_explained[0]*100:.1f}%)')
    plt.ylabel(f'Principal Component 2 ({var_explained[1]*100:.1f}%)')

    plt.tight_layout()
    plt.savefig(plot_path)
    print(f"PCA Plot saved to {plot_path}")

    if interactive:
        fig = px.scatter(pca_df, 
                        x='PC1', 
                        y='PC2', 
                        color='Group', 
                        hover_data=['Species', 'File'], 
                        color_discrete_map={'Endosymbionts': '#FC8D62', 'Free-Living Relatives': '#66C2A5'},
                        title='Interactive PCA of Genomic Features')
        
        fig.write_html(interactive_path)
        print(f"Interactive PCA Plot saved to {interactive_path}")

    loadings_df = pd.DataFrame(
        pca.components_.T,
        columns=['PC1_Loading', 'PC2_Loading'],
        index=feature_columns
    )
    loadings_df = loadings_df.sort_values('PC1_Loading', key=abs, ascending=False)
    loadings_df.to_csv(loadings_path)
    print(f"PCA loadings saved to {loadings_path}")

    return X, df[y], scaler

def run_random_forest(X, y_encoded, groups, le, outpath, n_splits=5, suffix=''):
    rf = RandomForestClassifier(
        n_estimators=500,
        max_depth=None,
        class_weight="balanced",
        random_state=28,
        n_jobs=-1,
    )

    cv = StratifiedGroupKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=28)

    y_pred = cross_val_predict(rf, X, y_encoded, groups=groups, cv=cv)
    y_prob = cross_val_predict(rf, X, y_encoded, groups=groups, cv=cv, method="predict_proba")

    acc = accuracy_score(y_encoded, y_pred)
    class_names = le.classes_
    print(f"  RF CV accuracy: {round(acc, 3)}")

    report = classification_report(y_encoded, y_pred, target_names=class_names)
    print('Classification Report:\n', report)
    with open(os.path.join(outpath, f'rf_classification_report{suffix}.txt'), "w") as f:
        f.write(report)

    # Confusion matrix
    cm = confusion_matrix(y_encoded, y_pred)
    fig, ax = plt.subplots(figsize=(6, 5))
    plt.rcParams.update({'font.size': 13})
    ConfusionMatrixDisplay(cm, display_labels=class_names).plot(ax=ax, cmap="Blues")
    ax.set_title(f"RF Confusion Matrix (CV Accuracy: {round(acc, 2)})")
    ax.set_xlabel('Predicted Label', fontweight='bold')
    ax.set_ylabel('True Label', fontweight='bold')
    plt.tight_layout()
    fig.savefig(os.path.join(outpath, f'rf_confusion_matrix{suffix}.pdf'))
    plt.close(fig)

    prob_df = pd.DataFrame({
        'File': X.index.get_level_values('File'),
        'Species': groups,
        'True_Label': le.inverse_transform(y_encoded),
        'Predicted_Label': le.inverse_transform(y_pred),
        'Endosymb_Probability': y_prob[:, 0]
    })
    prob_path = os.path.join(outpath, f'rf_prediction_probabilities{suffix}.csv')
    prob_df.to_csv(prob_path, index=False)
    print(f"Prediction probabilities saved to {prob_path}")

    # Flag high-confidence misclassifications
    prob_df['Confidence'] = prob_df['Endosymb_Probability'].apply(lambda p: max(p, 1 - p))
    errors_df = prob_df[prob_df['True_Label'] != prob_df['Predicted_Label']].copy()
    errors_df = errors_df.sort_values('Confidence', ascending=False)

    error_path = os.path.join(outpath, f'rf_high_confidence_errors{suffix}.csv')
    high_conf_errors = errors_df[errors_df['Confidence'] >= confidence_threshold]
    high_conf_errors.to_csv(error_path, index=False)

    print(f"\nHigh-Confidence Misclassifications (>= {confidence_threshold})")
    if high_conf_errors.empty:
        print("  None found.")
    else:
        print(high_conf_errors[['Species', 'File', 'True_Label', 'Predicted_Label', 'Endosymb_Probability', 'Confidence']].to_string(index=False))
    print(f"High-confidence misclassifications saved to {error_path}")

    # Fit on full data for importance
    rf.fit(X, y_encoded)
    # Gini importances
    label_map = {
        'Delta_GC2_4': '\u0394GC2-4',
        'GC4':         'GC4 Content',
        'AV_Bias':     'Ala/Val Bias',
        'Rest_Bias':   'Rest of AA Bias',
        'Transposase_Per_Gene': 'Transposase Density',
        'Mean_IGS_Size': 'Mean IGS Size',
    }
    imp = pd.DataFrame({
        'feature': [label_map.get(f, f) for f in feature_columns],
        'importance': rf.feature_importances_,
    }).sort_values('importance', ascending=True)

    fig, ax = plt.subplots(figsize=(8, 5))
    plt.rcParams.update({'font.size': 14})
    ax.barh(imp['feature'], imp['importance'])
    ax.set_xlabel('Gini Importance (Model Contribution)', fontweight='bold')
    ax.set_title('Feature Importance (Gini)', fontweight='bold')
    plt.tight_layout()
    fig.savefig(os.path.join(outpath, f'rf_feature_importance{suffix}.pdf'))
    plt.close(fig)

    # Permutation importances
    perm = permutation_importance(rf, X, y_encoded, n_repeats=30,
                                  random_state=28, n_jobs=-1)
    perm_df = pd.DataFrame({
        'feature': feature_columns,
        'importance_mean': perm.importances_mean,
        'importance_std': perm.importances_std,
    }).sort_values('importance_mean', ascending=False)
    perm_df.to_csv(os.path.join(outpath, f'rf_permutation_importance{suffix}.csv'), index=False)

    print('Random Forest analysis complete. Feature importance and permutation importance saved.')

    return rf

def rf_reduced(X, y_encoded, groups, le, outpath, n_splits=5):
    """Exclude recent endosymbionts and re-run the same model."""
    mask = ~groups.apply(is_recent_endosymb).values
    removed = sorted(groups[~mask].unique())
    print(f"\n  Excluding {(~mask).sum()} genomes from {len(removed)} recent endosymbiont species: {removed}")
    return run_random_forest(X[mask], y_encoded[mask], groups[mask], le, outpath, n_splits, suffix='_reduced')


if __name__ == "__main__":
    path = 'endosymb+relatives'
    input_path = os.path.join(path, feature_dir, features_df)
    outpath = os.path.join(path, feature_dir)

    plot_correlation(input_path, outpath)

    X, y, scaler = run_pca(input_path, outpath)

    plot_feature_distributions(X, y, outpath)

    le = LabelEncoder()
    y_encoded = le.fit_transform(y['Group'])
    groups = y['Species']

    rf_model = run_random_forest(X, y_encoded, groups, le, outpath)
    rf_reduced_model = rf_reduced(X, y_encoded, groups, le, outpath)

    # Serialize the model
    model_path = os.path.join('files', 'rf_model.joblib')
    os.makedirs('files', exist_ok=True)
    dump((rf_model, scaler, le), model_path)

    # Use on new genomes
    infer_path = os.path.join('ncbi_query', feature_dir, features_df)
    if os.path.exists(infer_path):
        rf_infer, scaler_infer, le_infer = load(model_path)
        unknown_df = pd.read_csv(infer_path)
        X_unk = unknown_df[feature_columns]
        X_unk_scaled = scaler_infer.transform(X_unk)
        probs = rf_infer.predict_proba(X_unk_scaled)
        preds = le_infer.inverse_transform(rf_infer.predict(X_unk_scaled))
        endosymb_idx = list(le_infer.classes_).index('endosymb_only')
        results_df = pd.DataFrame({
            'Species': unknown_df['Species'],
            'File': unknown_df['File'],
            'Predicted_Label': preds,
            'Endosymb_Probability': probs[:, endosymb_idx],
            'Confidence': [max(p, 1 - p) for p in probs[:, endosymb_idx]],
            'Low_Confidence': [max(p, 1 - p) < 0.6 for p in probs[:, endosymb_idx]],
        })
        results_path = os.path.join('ncbi_query', feature_dir, 'rf_predictions.csv')
        results_df.to_csv(results_path, index=False)
        print(f"Predictions for new genomes saved to {results_path}")