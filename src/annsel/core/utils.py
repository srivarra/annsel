from collections.abc import Iterable
from typing import Literal

import anndata as ad
import pandas as pd
from attr import dataclass
from narwhals.typing import IntoExpr
from scipy import sparse

from annsel.core import ObsExpr, ObsNamesExpr, VarExpr, VarNamesExpr, XExpr


def _extract_X(
    adata: ad.AnnData,
    layer: str | None = None,
    keep_sparse: bool = True,
    sparse_method: Literal["csr", "csc"] | None = None,
) -> pd.DataFrame:
    match (keep_sparse, sparse.issparse(adata.layers[layer]), sparse_method):
        case (True, True, _):
            _X = adata.layers[layer]
            return pd.DataFrame.sparse.from_spmatrix(data=_X, columns=adata.var_names, index=adata.obs_names)
        case (True, False, sparse_method):
            # Conver to sparse csr
            if sparse_method == "csr":
                _X = sparse.csr_matrix(adata.layers[layer])
            elif sparse_method == "csc":
                _X = sparse.csc_matrix(adata.layers[layer])
            else:
                raise ValueError(f"Invalid sparse_method: {sparse_method}")
            return pd.DataFrame.sparse.from_spmatrix(data=_X, columns=adata.var_names, index=adata.obs_names)
        case _:
            return adata.to_df(layer=layer)


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


# @nw.narwhalify
# def _compare_indices(first: IntoDataFrameT, second: IntoDataFrameT) -> IntoSeriesT:
#     return first & second


# def _fold_compare(*indices: pd.DataFrame) -> pd.Series:
#     start = indices[0]
#     for i in indices[1:]:
#         start = _compare_indices(start, i)
#     return start
