import os
import pandas as pd
import pymc as pm
import arviz as az
from bayesian_model_mvn import build_mvn_model

df = pd.read_csv('endosymb+relatives/feature_files/processed/processed_bayesian_data.csv')

feature_columns = ['GC4','AV_Bias','Rest_Bias',
                   'Transposase_Per_Gene','Mean_IGS_Size']

os.makedirs('traces', exist_ok=True)

# Rebuild indices (no subsampling)
df['clade_idx'] = df['Species'].astype('category').cat.codes
df['endo_idx'] = -1
mask = df['is_endosymbiont'] == 1
df.loc[mask, 'endo_idx'] = range(mask.sum())

model = build_mvn_model(df, feature_columns)
with model:
    trace = pm.sample(10000, tune=2000, chains=4,
                      target_accept=0.95, cores=1,
                      random_seed=42)

print(az.summary(trace, var_names=['mu_global','beta_decay',
                                    'sigma_clade']))
print(f"\nDivergences: {trace.sample_stats.diverging.sum().values}")

az.to_netcdf(trace, 'traces/full_mvn_trace.nc')
print("Trace saved.")