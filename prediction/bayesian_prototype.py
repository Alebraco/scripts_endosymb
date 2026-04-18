import pandas as pd
import pymc as pm
import arviz as az
import sys, os
from bayesian_model import build_prototype

df = pd.read_csv('endosymb+relatives/feature_files/processed/processed_bayesian_data.csv')

endo = df[df['is_endosymbiont'] == 1].sample(frac=125, random_state=42)
free = df[df['is_endosymbiont'] == 0].sample(frac=375, random_state=42)
sub = pd.concat([endo, free]).reset_index(drop=True)

sub['clade_idx'] = sub['Species'].astype('category').cat.codes
sub['endo_idx'] = -1
mask = sub['is_endosymbiont'] == 1
sub.loc[mask, 'endo_idx'] = range(mask.sum())

model = build_prototype(sub)
with model:
    trace = pm.sample(draws = 1000, tune=1000, target_accept=0.9,
                      chains=2, random_seed=42)

print(az.summary(trace, var_names=['mu_global', 'beta_decay', 'sigma_clade', 'sigma_obs']))
az.to_netcdf(trace, 'bayesian_prototype_trace.nc')