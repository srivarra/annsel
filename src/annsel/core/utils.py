from collections.abc import Iterable
from dataclasses import dataclass
from functools import reduce
from operator import and_
from typing import Literal

import anndata as ad
import pandas as pd
from narwhals.typing import IntoExpr
from scipy import sparse

from annsel.core import ObsExpr, ObsNamesExpr, VarExpr, VarNamesExpr, XExpr


def _extract_X(
    adata: ad.AnnData,
    layer: str | None = None,
    keep_sparse: bool = True,
    sparse_method: Literal["csr", "csc"] | None = None,
) -> pd.DataFrame:
    X = adata.layers[layer] if layer else adata.X
    match (keep_sparse, sparse.issparse(X), sparse_method):
        case (True, True, None):
            return pd.DataFrame.sparse.from_spmatrix(data=X, columns=adata.var_names, index=adata.obs_names)
        case (True, _, "csr"):
            return pd.DataFrame.sparse.from_spmatrix(
                data=sparse.csr_matrix(X), columns=adata.var_names, index=adata.obs_names
            )
        case (True, _, "csc"):
            return pd.DataFrame.sparse.from_spmatrix(
                data=sparse.csc_matrix(X), columns=adata.var_names, index=adata.obs_names
            )
        case (False, _, _):
            return pd.DataFrame(data=X.toarray(), columns=adata.var_names, index=adata.obs_names)
        case _:
            return pd.DataFrame(data=X, columns=adata.var_names, index=adata.obs_names)


@dataclass
class GroupedPredicates:
    """Grouped predicates."""

    obs: list[IntoExpr]
    var: list[IntoExpr]
    x: list[IntoExpr]
    obs_names: list[IntoExpr]
    var_names: list[IntoExpr]


def _map_predicates(
    *predicates: Iterable[IntoExpr | Iterable[IntoExpr] | list[bool]],
) -> GroupedPredicates:
    _obs = []
    _var = []
    _obs_names = []
    _var_names = []
    _x = []
    for predicate in predicates:
        match p := predicate[0] if isinstance(predicate, tuple) else predicate:
            case ObsExpr():
                _obs.append(p)
            case VarExpr():
                _var.append(p)
            case XExpr():
                _x.append(p)
            case ObsNamesExpr():
                _obs_names.append(p)
            case VarNamesExpr():
                _var_names.append(p)
            case _:
                raise ValueError(f"Invalid predicate: {predicate}")
    return GroupedPredicates(
        _obs,
        _var,
        _x,
        _obs_names,
        _var_names,
    )


def _get_final_indices(
    obj_names: pd.Index,
    *idx: pd.Index,
) -> pd.Index:
    return obj_names.intersection(pd.Index(list(reduce(and_, map(set, *idx)))))


def _handle_sparse_method(adata: ad.AnnData, sparse_method: Literal["csr", "csc"] | None) -> ad.AnnData:
    if sparse_method == "csc":
        _X = sparse.csc_matrix(adata.X)
    if sparse_method == "csr":
        _X = sparse.csr_matrix(adata.X)
    elif sparse_method is None:
        _X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
    return ad.AnnData(
        X=_X,
        obs=adata.obs,
        var=adata.var,
        obsm=adata.obsm,
        varm=adata.varm,
        obsp=adata.obsp,
        varp=adata.varp,
        layers=adata.layers,
        raw=adata.raw,
        uns=adata.uns,
    )
