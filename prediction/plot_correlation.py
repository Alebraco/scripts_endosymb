import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import plotly.express as px
import numpy as np
import os

feature_dir = 'feature_files'
features_df = 'combined_features.csv'

def plot_correlation(csv_path, outpath):
    plot_path = os.path.join(outpath, 'feature_correlation.pdf')
    df = pd.read_csv(csv_path)

    features = df.drop(columns=['Group', 'Species', 'File'])

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


def run_pca(csv_path, outpath):
    plot_path = os.path.join(outpath, 'pca_plot.pdf')
    interactive_path = os.path.join(outpath, 'pca_plot.html')
    df = pd.read_csv(csv_path)
    
    df['Group'] = df['Group'].replace({
        'endosymb_only': 'Endosymbionts',
        'relatives_only': 'Free-Living Relatives'
    })
    
    features = ['Delta_GC2_4', 'GC4', 'AV_Bias', 'Rest_Bias', 'Transposase_Per_Gene', 'Mean_IGS_Size']
    
    df_ml = df.set_index(['Group', 'Species', 'File'])
    X = df_ml[features]
    
    # Standardize data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(X_scaled)
    
    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'], index=df_ml.index)
    pca_df = pca_df.reset_index() 
    
    var_explained = pca.explained_variance_ratio_
    print(f"PC1 explains {var_explained[0]*100:.1f}% of the variance")
    print(f"PC2 explains {var_explained[1]*100:.1f}% of the variance")

    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        x='PC1', y='PC2',
        hue='Group',
        palette={'Endosymbionts': '#FC8D62', 'Free-Living Relatives': '#66C2A5'},
        data=pca_df,
        s=50, alpha=0.8, edgecolor='k'
    )

    plt.title('PCA of Genomic Features', fontsize=16, fontweight='bold')
    plt.legend(title='Group')
    plt.tight_layout()
    plt.savefig(plot_path)
    print(f"PCA Plot saved to {plot_path}")

    fig = px.scatter(pca_df, 
                     x='PC1', 
                     y='PC2', 
                     color='Group', 
                     hover_data=['Species', 'File'], 
                     color_discrete_map={'Endosymbionts': '#FC8D62', 'Free-Living Relatives': '#66C2A5'},
                     title='Interactive PCA of Genomic Features')
    
    fig.write_html(interactive_path)
    print(f"Interactive PCA Plot saved to {interactive_path}")

    return X, df['Group']

if __name__ == "__main__":
    path = 'endosymb+relatives'
    input_path = os.path.join(path, feature_dir, features_df)
    outpath = os.path.join(path, feature_dir)

    plot_correlation(input_path, outpath)
    X, y = run_pca(input_path, outpath)