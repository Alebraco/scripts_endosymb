import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import joblib
import os

df = pd.read_csv('endosymb+relatives/feature_files/combined_features.csv')

feature_columns = ['GC4', 'AV_Bias','Rest_Bias',
                   'Transposase_Per_Gene','Mean_IGS_Size']
scaler = StandardScaler()
df[feature_columns] = scaler.fit_transform(df[feature_columns])

df['Symbiotic_Status'] = df['Group'].map({'Endosymbiont': 1, 'Free-living': 0})
species_cat = df['Species'].astype('category')
df['clade_idx'] = species_cat.cat.codes

endosymb_mask = df['Symbiotic_Status'] == 1
df['endosymb_idx'] = -1
df.loc[endosymb_mask, 'endosymb_idx'] = np.arange(endosymb_mask.sum())

os.makedirs('endosymb+relatives/feature_files/processed', exist_ok=True)
df.to_csv('endosymb+relatives/feature_files/processed/processed_bayesian_data.csv', index=False)
joblib.dump(scaler, 'endosymb+relatives/feature_files/processed/bayesian_scaler.joblib')

print(f'Genomes: {len(df)}')
print(f'Endosymbionts: {df["Symbiotic_Status"].sum()}')
print(f'Clades: {df["clade_idx"].nunique()}')
print(f'Features: {feature_columns}')