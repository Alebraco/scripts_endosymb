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
from joblib import dump
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
from utils import feature_columns, POSTER_RCPARAMS, files_dir


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
    label_map = {
        'Delta_GC2_4':          'ΔGC2–4',
        'GC4':                  'GC4 Content',
        'AV_Bias':              'Codon usage bias (A/V)',
        'Rest_Bias':            'Codon usage bias (P/T/G/L/R/S)',
        'Transposase_Per_Gene': 'Transposase / Gene',
        'Mean_IGS_Size':        'Mean IGS Size',
    }

    groups_order = ['Endosymbionts', 'Free-Living Relatives']
    colors = ['#FC8D62', '#66C2A5']
    group_vals = y['Group'].values

    def _draw(figsize, xtick_fs):
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()
        for ax, feat in zip(axes, feature_columns):
            data = [X[feat].values[group_vals == g] for g in groups_order]
            bp = ax.boxplot(
                data,
                positions=range(1, len(groups_order) + 1),
                patch_artist=True,
                widths=0.5,
                flierprops=dict(marker='.', markersize=3, alpha=0.4),
                medianprops=dict(color='black', linewidth=1.5),
            )
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            ax.set_xticks(range(1, len(groups_order) + 1))
            if xtick_fs is not None:
                ax.set_xticklabels(groups_order, fontsize=xtick_fs)
            else:
                ax.set_xticklabels(groups_order)
            ax.set_title(label_map.get(feat, feat), fontweight='bold')
            ax.set_xlabel('')
        plt.tight_layout()
        return fig

    fig = _draw(figsize=(15, 10), xtick_fs=9)
    fig.savefig(os.path.join(outpath, 'feature_distributions_grid.pdf'))
    plt.close(fig)

    with plt.rc_context(POSTER_RCPARAMS):
        fig = _draw(figsize=(18, 11), xtick_fs=None)
        fig.savefig(os.path.join(outpath, 'fig_B_feature_boxplots.png'))
        plt.close(fig)

    print("Feature distribution plots saved.")

def run_pca(csv_path, outpath, label_recent=False, interactive=False, no_decorations=False, poster=False):
    plot_path = os.path.join(outpath, 'pca_plot_notitle.pdf' if (no_decorations and not label_recent) else 'pca_plot.pdf')
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

    if not (no_decorations and not label_recent):
        plt.title('PCA of Genomic Features', fontsize=16, fontweight='bold')
    if no_decorations and not label_recent:
        legend = plt.gca().get_legend()
        if legend:
            legend.remove()
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

    if poster:
        cluster_csv = os.path.join(files_dir, 'endosymbiont_clusters.csv')
        if not os.path.exists(cluster_csv):
            print(f"  Skipping fig_F: {cluster_csv} not found.")
        else:
            clusters = (
                pd.read_csv(cluster_csv)[['File', 'bgmm_prob', 'stage']]
                .dropna()
                .sort_values('bgmm_prob', ascending=False)
                .drop_duplicates(subset='File')
                [['File', 'stage']]
            )
            poster_df = pca_df.merge(clusters, on='File', how='left')
            poster_df.loc[poster_df['Group'] == 'Free-Living Relatives', 'stage'] = 'free-living'
            poster_df['Category'] = poster_df['stage'].fillna('free-living').str.capitalize()

            stage_viridis = plt.cm.viridis(np.linspace(0.2, 0.9, 4))
            palette = {
                'Recent':       stage_viridis[0],
                'Transitional': stage_viridis[1],
                'Reduced':      stage_viridis[2],
                'Ancient':      stage_viridis[3],
                'Free-living':  '#95A5A6',
            }
            draw_order = ['Free-living', 'Recent', 'Transitional', 'Reduced', 'Ancient']
            alpha_map = {'Free-living': 0.5}

            with plt.rc_context(POSTER_RCPARAMS):
                fig, ax = plt.subplots(figsize=(10, 8))
                for cat in draw_order:
                    sub = poster_df[poster_df['Category'] == cat]
                    if sub.empty:
                        continue
                    ax.scatter(
                        sub['PC1'], sub['PC2'],
                        color=palette[cat], label=cat,
                        s=70, alpha=alpha_map.get(cat, 0.85),
                        edgecolor='k', linewidths=0.3,
                    )
                ax.set_xlabel(f'Principal Component 1 ({var_explained[0]*100:.1f}%)')
                ax.set_ylabel(f'Principal Component 2 ({var_explained[1]*100:.1f}%)')
                ax.legend(loc='best')
                plt.tight_layout()
                fig.savefig(os.path.join(outpath, 'fig_F_pca_unified.png'))
                plt.close(fig)

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
    ConfusionMatrixDisplay(cm, display_labels=class_names).plot(ax=ax, cmap="Blues")
    ax.set_title(f"Confusion Matrix (CV Accuracy: {round(acc, 2)})", fontweight='bold')
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
        'Transposase_Per_Gene': 'Transposase\nDensity',
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

    # Poster composite: confusion matrix + feature importance side-by-side
    if suffix == '':
        with plt.rc_context(POSTER_RCPARAMS):
            fig, (ax_cm, ax_fi) = plt.subplots(1, 2, figsize=(16, 6))
            disp = ConfusionMatrixDisplay(cm, display_labels=class_names)
            disp.plot(ax=ax_cm, cmap='Blues', colorbar=False)
            for t in disp.text_.ravel():
                t.set_fontsize(26)
            ax_cm.set_title(f'Confusion Matrix (CV Accuracy: {round(acc, 2)})',
                            fontweight='bold')
            ax_cm.set_xlabel('Predicted Label', fontweight='bold')
            ax_cm.set_ylabel('True Label', fontweight='bold')

            ax_fi.barh(imp['feature'], imp['importance'])
            ax_fi.set_xlabel('Gini Importance (Model Contribution)',
                             fontweight='bold')
            ax_fi.set_title('Feature Importance (Gini)', fontweight='bold')

            plt.tight_layout()
            fig.savefig(os.path.join(outpath, 'fig_C_rf_classifier.png'))
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

    X, y, scaler = run_pca(input_path, outpath, poster=True)
    run_pca(input_path, outpath, no_decorations=True)

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