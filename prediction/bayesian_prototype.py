import pandas as pd
import pymc as pm
import arviz as az
from bayesian_model import build_prototype
import matplotlib.pyplot as plt

df = pd.read_csv('endosymb+relatives/feature_files/processed/processed_bayesian_data.csv')

endo = df[df['is_endosymbiont'] == 1].sample(n=125, random_state=42)
free = df[df['is_endosymbiont'] == 0].sample(n=375, random_state=42)
sub = pd.concat([endo, free]).reset_index(drop=True)

sub['clade_idx'] = sub['Species'].astype('category').cat.codes
sub['endo_idx'] = -1
mask = sub['is_endosymbiont'] == 1
sub.loc[mask, 'endo_idx'] = range(mask.sum())

model = build_prototype(sub)
with model:
    trace = pm.sample(draws=1000, tune=1000, target_accept=0.9,
                      chains=2, random_seed=42)

print(az.summary(trace, var_names=['mu_global', 'beta_decay', 'sigma_clade', 'sigma_obs']))
az.to_netcdf(trace, 'bayesian_prototype_trace.nc')

az.plot_trace(trace, var_names=['mu_global', 'beta_decay', 'sigma_clade', 'sigma_obs'])
plt.tight_layout()
plt.savefig('bayesian_prototype_trace.png')

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

cols = ['Species', 'File', 'Group', 'is_endosymbiont', 
        'theta_mean', 'theta_sd', 'GC4']
result[cols].to_csv('theta_estimates_prototype.csv', index=False)

endo_only = result[result['is_endosymbiont'] == 1]
print(endo_only[cols].tail(40))

print(result[cols].head(20))
print("---")
print(result[cols].tail(20))