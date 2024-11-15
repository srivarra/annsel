from collections.abc import Iterable
from typing import Literal

import anndata as ad
import narwhals as nw
import pandas as pd
from narwhals.typing import IntoDataFrame, IntoExpr

from annsel.core.models import NarwhalsMethod, predicate_guard
from annsel.core.typing import XIndicies
from annsel.core.utils import _extract_X, _get_final_indices, _map_predicates


@nw.narwhalify
def _filter_df(df: IntoDataFrame, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]) -> IntoDataFrame:
    return df.filter(*predicates)


class FilterAnnData(NarwhalsMethod):
    def __init__(self, adata: ad.AnnData, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]):
        self._adata = adata
        self._predicates = _map_predicates(*predicates)

    def _apply_predicates(self, df: pd.DataFrame, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]) -> pd.Index:
        return _filter_df(df, *predicates).index

    @predicate_guard("var")
    def _run_var_predicates(self) -> pd.Index:
        return self._apply_predicates(self._adata.var, *self._predicates.var)

    @predicate_guard("var_names")
    def _run_var_names_predicates(self) -> pd.Index:
        return self._apply_predicates(self._adata.var_names.to_frame(name="var_names"), *self._predicates.var_names)

    @predicate_guard("obs")
    def _run_obs_predicates(self) -> pd.Index:
        return self._apply_predicates(self._adata.obs, *self._predicates.obs)

    @predicate_guard("obs_names")
    def _run_obs_names_predicates(self) -> pd.Index:
        return self._apply_predicates(self._adata.obs_names.to_frame(name="obs_names"), *self._predicates.obs_names)

    @predicate_guard("x")
    def _run_x_predicates(
        self, layer: str | None = None, keep_sparse: bool = True, sparse_method: Literal["csr", "csc"] | None = None
    ) -> XIndicies:
        _X_df: pd.DataFrame = _extract_X(
            adata=self._adata, layer=layer, keep_sparse=keep_sparse, sparse_method=sparse_method
        )

        _X_df = _filter_df(_X_df, *self._predicates.x)
        return XIndicies(_X_df.index, _X_df.columns)

    def _finalize_indices_obs(self, *idx: pd.Index) -> pd.Index:
        return _get_final_indices(self._adata.obs_names, *idx)

    def _finalize_indices_var(self, *idx: pd.Index) -> pd.Index:
        return _get_final_indices(self._adata.var_names, *idx)

    def __call__(
        self,
        layer: str | None = None,
        keep_sparse: bool = True,
        sparse_method: Literal["csr", "csc"] | None = None,
    ) -> tuple[pd.Index, pd.Index]:
        obs_indices = []
        var_indices = []

        # Run predicates - guards will return None if no predicates exist
        if (var_idx := self._run_var_predicates()) is not None:
            var_indices.append(var_idx)
        if (var_names_idx := self._run_var_names_predicates()) is not None:
            var_indices.append(var_names_idx)
        if (obs_idx := self._run_obs_predicates()) is not None:
            obs_indices.append(obs_idx)
        if (obs_names_idx := self._run_obs_names_predicates()) is not None:
            obs_indices.append(obs_names_idx)
        if x_indices := self._run_x_predicates(layer=layer, keep_sparse=keep_sparse, sparse_method=sparse_method):
            obs_indices.append(x_indices.obs)
            var_indices.append(x_indices.var)

        # Finalize indices or use defaults
        final_obs_idx = self._finalize_indices_obs(obs_indices) if obs_indices else self._adata.obs_names
        final_var_idx = self._finalize_indices_var(var_indices) if var_indices else self._adata.var_names

        return final_obs_idx, final_var_idx
