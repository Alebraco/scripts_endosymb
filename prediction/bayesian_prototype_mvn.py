import pandas as pd
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt
from bayesian_model_mvn import build_mvn_model

df = pd.read_csv('endosymb+relatives/feature_files/processed/processed_bayesian_data.csv')

feature_columns = ['GC4','AV_Bias','Rest_Bias',
                   'Transposase_Per_Gene','Mean_IGS_Size']

# Subsample 500 genomes
endo = df[df['is_endosymbiont']==1].sample(125, random_state=42)
free = df[df['is_endosymbiont']==0].sample(375, random_state=42)
sub = pd.concat([free, endo]).reset_index(drop=True)

# Rebuild indices
sub['clade_idx'] = sub['Species'].astype('category').cat.codes
sub['endo_idx'] = -1
mask = sub['is_endosymbiont']==1
sub.loc[mask, 'endo_idx'] = range(mask.sum())

# Build and sample
model = build_mvn_model(sub, feature_columns)
with model:
    trace = pm.sample(1000, tune=1000, chains=2,
                      target_accept=0.9, random_seed=42)

# Diagnostics
print(az.summary(trace, var_names=['mu_global','beta_decay',
                                    'sigma_clade']))
print(f"\nDivergences: {trace.sample_stats.diverging.sum().values}")

az.to_netcdf(trace, 'mvn_prototype_trace.nc')

# Extract theta summaries
theta_endo_summary = az.summary(trace, var_names=['theta_endo'])
theta_free_summary = az.summary(trace, var_names=['theta_free'])

endo_rows = sub[sub['is_endosymbiont'] == 1].copy()
endo_rows['theta_mean'] = theta_endo_summary['mean'].values
endo_rows['theta_sd'] = theta_endo_summary['sd'].values

free_rows = sub[sub['is_endosymbiont'] == 0].copy()
free_rows['theta_mean'] = theta_free_summary['mean'].values
free_rows['theta_sd'] = theta_free_summary['sd'].values

result = pd.concat([endo_rows, free_rows])
result = result.sort_values('theta_mean', ascending=False)

endo_only = result[result['is_endosymbiont'] == 1]
print("TOP 10:")
print(endo_only[['Species','theta_mean','theta_sd']].head(10))
print("\nBOTTOM 10:")
print(endo_only[['Species','theta_mean','theta_sd']].tail(10))

result.to_csv('theta_estimates_mvn_prototype.csv', index=False)

# Trace plots
# Global parameters
az.plot_trace(trace, var_names=['mu_global','beta_decay',
                                 'sigma_clade'])
plt.tight_layout()
plt.savefig('trace_plots_mvn_global.png')
plt.close()

# Sample of theta values (not all 125)
az.plot_trace(trace, var_names=['theta_endo'], 
              coords={'theta_endo_dim_0': range(10)})
plt.tight_layout()
plt.savefig('trace_plots_mvn_theta_sample.png')
plt.close()

print("Trace plots saved.")