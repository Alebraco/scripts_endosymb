import pymc as pm
import numpy as np
import pytensor.tensor as pt

def build_mvn_model(df, feature_columns):
    """
    5-feature multivariate Bayesian hierarchical model.
    """
    
    # Data
    y = df[feature_columns].values
    clade_idx = df['clade_idx'].values
    is_endo = df['is_endosymbiont'].values
    
    n_features = len(feature_columns)
    n_clades = df['clade_idx'].nunique()
    n_endo = int(is_endo.sum())
    n_free = len(df) - n_endo
    endo_positions = np.where(is_endo == 1)[0]
    free_positions = np.where(is_endo == 0)[0]
    
    with pm.Model() as model:
        # Global parameters
        mu_global = pm.Normal('mu_global', mu=0, sigma=1, 
                              shape=n_features)
        beta_decay = pm.Normal('beta_decay', mu=0, sigma=1, 
                               shape=n_features)
        
        # Species-level random effects
        sigma_clade = pm.HalfNormal('sigma_clade', sigma=1, 
                                     shape=n_features)
        mu_clade = pm.Normal('mu_clade', mu=mu_global, 
                             sigma=sigma_clade, 
                             shape=(n_clades, n_features))
        
        # Transition parameters
        theta_endo = pm.HalfNormal('theta_endo', sigma=1, 
                                    shape=n_endo)
        theta_free = pm.HalfNormal('theta_free', sigma=0.01, 
                                    shape=n_free)
        
        # Build full theta
        theta = pt.zeros(len(y))
        theta = pt.set_subtensor(theta[free_positions], theta_free)
        theta = pt.set_subtensor(theta[endo_positions], theta_endo)
        
        # Expected value
        # mu_clade[clade_idx] is (N, 5)
        # theta[:, None] is (N, 1) * beta_decay[None, :] is (1, 5) -> (N, 5)
        mu = mu_clade[clade_idx] + theta[:, None] * beta_decay[None, :]
        
        # Covariance: LKJ prior
        chol, corr, stds = pm.LKJCholeskyCov(
            'packed_L',
            n=n_features,
            eta=2.0,
            sd_dist=pm.HalfNormal.dist(sigma=1),
            compute_corr=True
        )
        
        # Likelihood
        obs = pm.MvNormal('obs', mu=mu, chol=chol, observed=y)
    
    return model