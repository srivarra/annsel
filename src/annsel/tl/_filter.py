from typing import Literal

import anndata as ad
import narwhals as nw
import pandas as pd
from narwhals.typing import IntoDataFrameT

from annsel.core.typing import Predicates
from annsel.core.utils import _construct_adata_from_indices, _get_final_indices


@nw.narwhalify
def _filter_df(df: IntoDataFrameT, *predicates: Predicates) -> IntoDataFrameT:
    return df.filter(*predicates)


def _filter_obs(adata: ad.AnnData, *predicates: Predicates) -> pd.Index:
    return _filter_df(adata.obs, *predicates).index


def _filter_var(adata: ad.AnnData, *predicates: Predicates) -> pd.Index:
    return _filter_df(adata.var, *predicates).index


def _filter_x(adata: ad.AnnData, *predicates: Predicates, layer: str | None = None) -> pd.Index:
    return _filter_df(adata.to_df(layer=layer), *predicates).index


def _filter_obs_names(adata: ad.AnnData, *predicates: Predicates) -> pd.Index:
    return _filter_df(adata.obs_names.to_frame(name="obs_names"), *predicates).index


def _filter_var_names(adata: ad.AnnData, *predicates: Predicates) -> pd.Index:
    return _filter_df(adata.var_names.to_frame(name="var_names"), *predicates).index


def _filter(
    adata: ad.AnnData,
    obs: Predicates | None = None,
    var: Predicates | None = None,
    x: Predicates | None = None,
    obs_names: Predicates | None = None,
    var_names: Predicates | None = None,
    layer: str | None = None,
    sparse: Literal["csr", "csc", True, False] | None = None,
) -> ad.AnnData:
    obs_names_idx = []
    var_names_idx = []

    if obs:
        obs_names_idx.append(_filter_obs(adata, obs))

    if obs_names:
        obs_names_idx.append(_filter_obs_names(adata, obs_names))

    if var:
        var_names_idx.append(_filter_var(adata, var))

    if var_names:
        var_names_idx.append(_filter_var_names(adata, var_names))

    if x:
        obs_names_idx.append(_filter_x(adata, x, layer=layer))

    final_obs_idx = adata.obs_names if not obs_names_idx else _get_final_indices(adata.obs_names, obs_names_idx)
    final_var_idx = adata.var_names if not var_names_idx else _get_final_indices(adata.var_names, var_names_idx)

    return _construct_adata_from_indices(adata, final_obs_idx, final_var_idx, sparse)
