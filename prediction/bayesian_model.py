import pymc as pm
import numpy as np

def build_prototype(df, feature='GC4'):
    """Single-feature hierarchical model for transition parameter estimation."""
    
    y = df[feature].values
    clade_idx = df['clade_idx'].values
    is_endo = df['is_endosymbiont'].values

    n_clades = df['clade_idx'].nunique()
    n_endo = int(is_endo.sum())
    n_free = len(df) - n_endo
    endo_positions = np.where(is_endo == 1)[0]
    free_positions = np.where(is_endo == 0)[0]
    
    with pm.Model() as model:
        # Global parameters
        mu_global = pm.Normal('mu_global', mu=0, sigma=1)
        beta_decay = pm.Normal('beta_decay', mu=0, sigma=1)
        
        # Species-level random effects
        sigma_clade = pm.HalfNormal('sigma_clade', sigma=1)
        mu_clade = pm.Normal('mu_clade', mu=mu_global, 
                             sigma=sigma_clade, shape=n_clades)
        
        # Endosymbiont-specific decay signal
        theta_endo = pm.HalfNormal('theta_endo', sigma=1, 
                                    shape=n_endo)
        theta_free = pm.HalfNormal('theta_free', mu=0, sigma=0.01,
                                  shape=n_free)
        
        # Build full theta
        theta = pm.math.zeros_like(y, dtype='floatX')
        theta = pm.math.set_subtensor(
            theta[endo_positions], theta_endo)
        theta = pm.math.set_subtensor(
            theta[free_positions], theta_free)
        
        # Expected value
        mu = mu_clade[clade_idx] + theta * beta_decay
        
        # Likelihood
        sigma_obs = pm.HalfNormal('sigma_obs', sigma=1)
        obs = pm.Normal('obs', mu=mu, sigma=sigma_obs, 
                        observed=y)
    
    return model