import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import os

def plot_correlation(path):
    feature_dir = 'feature_files'
    input_path = os.path.join(path, feature_dir, 'combined_features.csv')
    plot_path = os.path.join(path, feature_dir, 'feature_correlation.pdf')
    
    df = pd.read_csv(input_path)
    features = df.drop(columns=['Group', 'Species', 'File'])

    corr = features.corr()

    plt.figure(figsize=(10, 8))
    sns.heatmap(corr, annot=True, cmap='coolwarm', fmt=".2f", linewidths=0.5)

    plt.title("Feature Correlation Matrix")
    plt.tight_layout()

    plt.savefig(plot_path)
    plt.close()
    
    upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))
    to_drop = [column for column in upper.columns if any(upper[column] > 0.90)]
    if to_drop:
        print(f"Highly correlated features (>0.9): {to_drop}")

if __name__ == "__main__":
    path = 'endosymb+relatives'
    plot_correlation(path)