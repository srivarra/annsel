import itertools
from collections.abc import Iterable, Iterator

import anndata as ad
import narwhals as nw
import pandas as pd
from narwhals.group_by import GroupBy as NwGroupBy
from narwhals.typing import IntoDataFrame, IntoExpr
from pandas import Index

from annsel.core.models import NarwhalsMethod, predicate_guard
from annsel.core.utils import _map_predicates


@nw.narwhalify
def _groupby_df(df: IntoDataFrame, *exprs: IntoExpr) -> NwGroupBy:
    return df.group_by(*exprs)


class GroupByAnnData(NarwhalsMethod):
    def __init__(self, adata: ad.AnnData, *exprs: IntoExpr | Iterable[IntoExpr] | list[bool]):
        self._adata = adata
        self._predicates = _map_predicates(*exprs)

    def _apply_predicates(self, df: pd.DataFrame, *exprs: IntoExpr | Iterable[IntoExpr] | list[bool]) -> NwGroupBy:
        return _groupby_df(df, *exprs)

    @predicate_guard("obs")
    def _run_obs_predicates(self) -> NwGroupBy:
        return self._apply_predicates(self._adata.obs, *self._predicates.obs)

    @predicate_guard("var")
    def _run_var_predicates(self) -> NwGroupBy:
        return self._apply_predicates(self._adata.var, *self._predicates.var)

    @predicate_guard("x")
    def _run_x_predicates(self) -> NwGroupBy:
        raise ValueError("X is not supported for `.groupby`")

    @predicate_guard("var_names")
    def _run_var_names_predicates(self) -> NwGroupBy:
        raise ValueError("var_names are not supported for `.groupby`")

    @predicate_guard("obs_names")
    def _run_obs_names_predicates(self) -> NwGroupBy:
        raise ValueError("obs_names are not supported for `.groupby`")

    def _finalize_indices_obs(self, *idx: Index) -> Index:
        pass

    def _finalize_indices_var(self, *idx: Index) -> Index:
        pass

    def __call__(self, return_group_names: bool = False) -> Iterator[ad.AnnData]:
        obs_gb_list = []
        var_gb_list = []

        if (var_gb := self._run_var_predicates()) is not None:
            var_gb_list.append(var_gb)
        if (obs_gb := self._run_obs_predicates()) is not None:
            obs_gb_list.append(obs_gb)
        if (self._run_x_predicates()) is not None:
            raise ValueError("X is not supported for `.groupby`")
        if (self._run_var_names_predicates()) is not None:
            raise ValueError("var_names are not supported for `.groupby`")
        if (self._run_obs_names_predicates()) is not None:
            raise ValueError("obs_names are not supported for `.groupby`")

        # Get the first (and only) group by object from each list
        obs_gb = obs_gb_list[0] if obs_gb_list else None
        var_gb = var_gb_list[0] if var_gb_list else None

        if obs_gb and not var_gb:
            for group_key, group_df in obs_gb:
                group_data = self._adata[group_df.to_native().index, :]
                if return_group_names:
                    yield (group_key, group_data)
                else:
                    yield group_data
        elif var_gb and not obs_gb:
            for group_key, group_df in var_gb:
                group_data = self._adata[:, group_df.to_native().index]
                if return_group_names:
                    yield (group_key, group_data)
                else:
                    yield group_data
        elif obs_gb and var_gb:
            for (obs_key, obs_df), (var_key, var_df) in itertools.product(obs_gb, var_gb):
                group_data = self._adata[obs_df.to_native().index, var_df.to_native().index]
                if return_group_names:
                    yield ((obs_key, var_key), group_data)
                else:
                    yield group_data
