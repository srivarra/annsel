from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

import anndata as ad
import narwhals as nw
import pandas as pd
from narwhals.typing import IntoDataFrame, IntoExpr

from annsel.core.utils import _extract_X


@nw.narwhalify
def _filter_df(df: IntoDataFrame, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]) -> IntoDataFrame:
    return df.filter(*predicates)


def _filter_adata_by_obs(adata: ad.AnnData, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]) -> pd.Index:
    return _filter_df(adata.obs, *predicates).index


def _filter_adata_by_var(adata: ad.AnnData, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]) -> pd.Index:
    return _filter_df(adata.var, *predicates).index


def _filter_adata_by_var_names(adata: ad.AnnData, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]) -> pd.Index:
    return _filter_df(adata.var_names.to_frame(name="var_names"), *predicates).index


def _filter_adata_by_obs_names(adata: ad.AnnData, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]) -> pd.Index:
    return _filter_df(adata.obs_names.to_frame(name="obs_names"), *predicates).index


@dataclass
class FilteredXIndicies:
    obs: pd.Index
    var: pd.Index


def _filter_adata_by_X(
    adata: ad.AnnData,
    *predicates: IntoExpr | Iterable[IntoExpr] | list[bool],
    layer=str | None,
    keep_sparse: bool = True,
    sparse_method: Literal["csr", "csc"] | None = None,
) -> FilteredXIndicies:
    _X_df: pd.DataFrame = _extract_X(adata, layer=layer, keep_sparse=keep_sparse, sparse_method=sparse_method)

    _X_df = _filter_df(_X_df, *predicates)
    return FilteredXIndicies(_X_df.index, _X_df.columns)
