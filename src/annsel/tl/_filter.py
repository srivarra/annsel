from collections.abc import Mapping

import anndata as ad
import narwhals as nw
import pandas as pd
from narwhals.typing import Frame

from annsel.core.typing import Predicates
from annsel.core.utils import _construct_adata_from_indices, _get_final_indices


@nw.narwhalify
def _filter_df(df: Frame, *predicates: Predicates) -> Frame:
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


def _filter_obsm(adata: ad.AnnData, key: str, *predicates: Predicates) -> pd.Index:
    return _filter_df(pd.DataFrame(adata.obsm[key], index=adata.obs_names), *predicates).index


def _filter_varm(adata: ad.AnnData, key: str, *predicates: Predicates) -> pd.Index:
    return _filter_df(pd.DataFrame(adata.varm[key], index=adata.var_names), *predicates).index


def _filter(
    adata: ad.AnnData,
    obs: Predicates | None = None,
    var: Predicates | None = None,
    x: Predicates | None = None,
    obs_names: Predicates | None = None,
    var_names: Predicates | None = None,
    obsm: Mapping[str, Predicates] | None = None,
    varm: Mapping[str, Predicates] | None = None,
    layer: str | None = None,
) -> ad.AnnData:
    obs_names_idx = []
    var_names_idx = []

    if obs is not None:
        obs_names_idx.append(_filter_obs(adata, obs))

    if obs_names is not None:
        obs_names_idx.append(_filter_obs_names(adata, obs_names))

    if var is not None:
        var_names_idx.append(_filter_var(adata, var))

    if var_names is not None:
        var_names_idx.append(_filter_var_names(adata, var_names))

    if x is not None:
        obs_names_idx.append(_filter_x(adata, x, layer=layer))

    if obsm is not None:
        for key, predicates in obsm.items():
            obs_names_idx.append(_filter_obsm(adata, key, predicates))

    if varm is not None:
        for key, predicates in varm.items():
            var_names_idx.append(_filter_varm(adata, key, predicates))

    final_obs_idx = adata.obs_names if not obs_names_idx else _get_final_indices(adata.obs_names, *obs_names_idx)
    final_var_idx = adata.var_names if not var_names_idx else _get_final_indices(adata.var_names, *var_names_idx)
    return _construct_adata_from_indices(adata, final_obs_idx, final_var_idx)
