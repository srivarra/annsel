from typing import Literal

import anndata as ad
import narwhals as nw
import pandas as pd
from narwhals.typing import IntoDataFrame

from annsel.core.typing import Predicates
from annsel.core.utils import _construct_adata_from_indices, _get_final_indices


@nw.narwhalify
def _select_df(df: IntoDataFrame, *predicates: Predicates) -> IntoDataFrame:
    return df.select(*predicates)


def _select_obs(adata: ad.AnnData, *predicates: Predicates) -> pd.Index:
    return _select_df(adata.obs, *predicates).columns


def _select_var(adata: ad.AnnData, *predicates: Predicates) -> pd.Index:
    return _select_df(adata.var, *predicates).columns


def _select_x(adata: ad.AnnData, *predicates: Predicates, layer: str | None = None) -> pd.Index:
    return _select_df(adata.to_df(layer=layer), *predicates).columns


def _select(
    adata: ad.AnnData,
    obs: Predicates | None = None,
    var: Predicates | None = None,
    x: Predicates | None = None,
    sparse: Literal["csr", "csc", True, False] | None = None,
) -> ad.AnnData:
    obs_columns = []
    var_columns = []
    var_names = []

    if obs:
        obs_columns.append(_select_obs(adata, obs))

    if var:
        var_columns.append(_select_var(adata, var))

    if x:
        var_names.append(_select_x(adata, x, layer=None))

    final_obs_cols = _get_final_indices(adata.obs.columns, obs_columns) if obs_columns else adata.obs.columns
    final_var_cols = _get_final_indices(adata.var.columns, var_columns) if var_columns else adata.var.columns
    final_var_names = _get_final_indices(adata.var_names, var_names) if var_names else adata.var_names

    _adata = _construct_adata_from_indices(
        adata, obs_idx=adata.obs_names, var_idx=final_var_names, sparse_method=sparse
    )

    _adata.obs = _adata.obs[final_obs_cols]
    _adata.var = _adata.var[final_var_cols]

    return _adata
