"""Filter operations using narwhals-anndata plugin backend.

This module provides the new implementation of filter operations
that leverage the narwhals-anndata plugin for better performance
and maintainability.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import anndata as ad

    from annsel.core.typing import Predicates


def _filter_plugin(
    adata: ad.AnnData,
    obs: Predicates | None = None,
    var: Predicates | None = None,
    x: Predicates | None = None,
    obs_names: Predicates | None = None,
    var_names: Predicates | None = None,
    layer: str | None = None,
) -> ad.AnnData:
    """Filter AnnData using the narwhals-anndata plugin.

    This is the new implementation that uses the narwhals-anndata plugin
    as the backend, providing better performance and ecosystem integration.

    Parameters
    ----------
    adata
        The AnnData object to filter.
    obs
        Predicates to filter observations by.
    var
        Predicates to filter variables by.
    x
        Predicates to filter by expression values in X matrix.
    obs_names
        Predicates to filter by observation names.
    var_names
        Predicates to filter by variable names.
    layer
        The layer to use for X filtering.

    Returns
    -------
    ad.AnnData
        The filtered AnnData object.
    """
    import narwhals as nw
    import pandas as pd
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    # Get plugin namespace with current narwhals version
    namespace = __narwhals_namespace__(Version.MAIN)
    result_adata = adata

    # Filter observations using plugin
    if obs is not None:
        nwdata = namespace.from_native(result_adata)
        result_adata = nwdata.obs.filter(obs)

    # Filter by observation names
    if obs_names is not None:
        obs_names_df = pd.DataFrame({"obs_names": result_adata.obs_names}, index=result_adata.obs_names)
        obs_names_nw = nw.from_native(obs_names_df)
        filtered_obs_names = obs_names_nw.filter(obs_names)
        filtered_indices = filtered_obs_names.to_native().index
        result_adata = result_adata[filtered_indices, :]

    # Filter variables using plugin
    if var is not None:
        nwdata = namespace.from_native(result_adata)
        result_adata = nwdata.var.filter(var)

    # Filter by variable names
    if var_names is not None:
        var_names_df = pd.DataFrame({"var_names": result_adata.var_names}, index=result_adata.var_names)
        var_names_nw = nw.from_native(var_names_df)
        filtered_var_names = var_names_nw.filter(var_names)
        filtered_indices = filtered_var_names.to_native().index
        result_adata = result_adata[:, filtered_indices]

    # Filter by X values (materialize to DataFrame)
    if x is not None:
        x_df = result_adata.to_df(layer=layer)
        x_nw = nw.from_native(x_df)
        filtered_x = x_nw.filter(x)
        filtered_indices = filtered_x.to_native().index
        result_adata = result_adata[filtered_indices, :]

    return result_adata
