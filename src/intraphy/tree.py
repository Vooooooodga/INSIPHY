"""Compatible tree API. Topology and legacy models have separate owners."""
from .topology import SpeciesTree

# Old public numeric helpers remain lazily importable for legacy clients.
_LEGACY_NAMES = frozenset({'discrete_likelihood_tables', 'ctmc_posteriors', '_transition_probability_approx', 'discrete_log_likelihood', 'transition_probability', 'chi_square_sf_df1', '_tip_observations', 'sample_ctmc_bridge', 'sankoff', '_quantile', 'fit_invariant_test', 'stochastic_map_summary', 'bootstrap_invariant_test', '_sample_weighted', 'transition_matrix', 'fit_foreground_rate_test', 'rate_matrix_array', 'fit_discrete_ctmc', 'transition_cost', '_transition_matrix_cached', 'emission_probability'})

def __getattr__(name):
    if name in _LEGACY_NAMES:
        from intraphy.experimental import tree_model
        return getattr(tree_model, name)
    raise AttributeError(name)
