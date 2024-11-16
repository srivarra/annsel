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
def _select_df(df: IntoDataFrame, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]) -> IntoDataFrame:
    return df.select(*predicates)


class SelectAnnData(NarwhalsMethod):
    def __init__(self, adata: ad.AnnData, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]):
        self._adata = adata
        self._predicates = _map_predicates(*predicates)

    def _apply_predicates(
        self,
        df: pd.DataFrame,
        *predicates: nw.Expr | str | nw.Series | Iterable[nw.Expr | str | nw.Series] | list[bool],
    ) -> pd.Index:
        return _select_df(df, *predicates).columns

    @predicate_guard("var")
    def _run_var_predicates(self) -> pd.Index:
        return self._apply_predicates(self._adata.var, *self._predicates.var)

    @predicate_guard("obs")
    def _run_obs_predicates(self) -> pd.Index:
        return self._apply_predicates(self._adata.obs, *self._predicates.obs)

    @predicate_guard("x")
    def _run_x_predicates(
        self, layer: str | None = None, keep_sparse: bool = True, sparse_method: Literal["csr", "csc"] | None = None
    ) -> XIndicies | None:
        _X_df: pd.DataFrame = _extract_X(self._adata, layer=layer, keep_sparse=keep_sparse, sparse_method=sparse_method)
        _X_cols = self._apply_predicates(_X_df, *self._predicates.x)
        return XIndicies(_X_df.index, _X_cols)

    @predicate_guard("var_names")
    def _run_obs_names_predicates(self) -> pd.Index:
        raise ValueError("var_names are not supported for `.select`, only `obs`, `var` and `x` are supported.")

    @predicate_guard("obs_names")
    def _run_obs_names_predicates(self) -> pd.Index:
        raise ValueError("obs_names are not supported for `.select`, only `obs`, `var` and `x` are supported.")

    def _finalize_indices_obs(self, *idx: pd.Index) -> pd.Index:
        return _get_final_indices(self._adata.obs.columns, *idx)

    def _finalize_indices_var(self, *idx: pd.Index) -> pd.Index:
        return _get_final_indices(self._adata.var.columns, *idx)

    def _finalize_var_names(self, *idx: pd.Index) -> pd.Index:
        return _get_final_indices(self._adata.var_names, *idx)

    @predicate_guard("var_names")
    def _run_var_names_predicates(self) -> pd.Index:
        raise ValueError(
            "var_names predicates are not supported for `.select`, only `obs`, `var` and `x` are supported."
        )

    @predicate_guard("obs_names")
    def _run_obx_names_predicates(self) -> pd.Index:
        raise ValueError(
            "var_names predicates are not supported for `.select`, only `obs`, `var` and `x` are supported."
        )

    def __call__(
        self,
        layer: str | None = None,
        keep_sparse: bool = True,
        sparse_method: Literal["csr", "csc"] | None = None,
    ) -> XIndicies:
        obs_columns = []
        var_columns = []
        var_names = []

        # Run predicates - guards will return None if no predicates exist
        if (var_col := self._run_var_predicates()) is not None:
            var_columns.append(var_col)
        if (obs_col := self._run_obs_predicates()) is not None:
            obs_columns.append(obs_col)
        if x_indices := self._run_x_predicates(layer=layer, keep_sparse=keep_sparse, sparse_method=sparse_method):
            var_names.append(x_indices.var)

        # Error out if the users pass these in
        self._run_var_names_predicates()
        self._run_obs_names_predicates()

        # Finalize indices or use defaults
        final_obs_cols = self._finalize_indices_obs(obs_columns) if obs_columns else self._adata.obs.columns
        final_var_cols = self._finalize_indices_var(var_columns) if var_columns else self._adata.var.columns
        final_var_names = self._finalize_var_names(var_names) if var_names else self._adata.var_names

        return final_obs_cols, final_var_cols, final_var_names
