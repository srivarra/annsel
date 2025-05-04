import itertools
from collections.abc import Callable, Generator
from dataclasses import dataclass
from functools import cached_property

import anndata as ad
import narwhals as nw
import pandas as pd
from more_itertools import first
from narwhals.typing import Frame

from annsel.core.typing import NwGroupBy, Predicates
from annsel.core.utils import _construct_adata_from_indices, _extract_names_from_expr, second


@nw.narwhalify
def _group_by_df(df: Frame, *names: str) -> NwGroupBy:
    return df.group_by(*names)


def _group_by_obs(adata: ad.AnnData, *predicates: Predicates):
    _gb = _group_by_df(adata.obs, *predicates)
    return (predicates, _gb)


def _group_by_var(adata: ad.AnnData, *predicates: Predicates):
    _gb = _group_by_df(adata.var, *predicates)
    return (predicates, _gb)


@dataclass
class GroupByAnndata:
    """Container for grouped AnnData objects and their metadata.

    Parameters
    ----------
    obs_values
        Values for the observation grouping columns in this group.
    var_values
        Values for the variable grouping columns in this group.
    adata
        The actual data subset.
    _obs_predicates
        Predicates used for observation grouping (internal).
    _var_predicates
        Predicates used for variable grouping (internal).

    Returns
    -------
    A container for grouped AnnData objects and their metadata.
    """

    # Values for those columns in this group
    obs_values: tuple[str, ...]
    var_values: tuple[str, ...]

    # The actual data subset
    adata: ad.AnnData

    # Predicates used for grouping (internal, not shown in repr)
    _obs_predicates: Predicates | None = None
    _var_predicates: Predicates | None = None

    @cached_property
    def _obs_column_names(self) -> tuple[str, ...]:
        """Cached property for observation column names."""
        if self._obs_predicates is None:
            return ()
        if isinstance(self._obs_predicates, str | nw.Expr):
            return _extract_names_from_expr(self._obs_predicates)
        else:
            return _extract_names_from_expr(*self._obs_predicates)

    @cached_property
    def _var_column_names(self) -> tuple[str, ...]:
        """Cached property for variable column names."""
        if self._var_predicates is None:
            return ()
        if isinstance(self._var_predicates, str | nw.Expr):
            return _extract_names_from_expr(self._var_predicates)
        else:
            return _extract_names_from_expr(*self._var_predicates)

    @cached_property
    def obs_dict(self):
        """Dictionary mapping observation column names to their values."""
        return dict(zip(self._obs_column_names, self.obs_values, strict=False))

    @cached_property
    def var_dict(self):
        """Dictionary mapping variable column names to their values."""
        var_cols = self._var_column_names
        return dict(zip(var_cols, self.var_values, strict=False))

    def _format_repr_branch(self, indent: str, branch_name: str, cols: tuple[str, ...], vals: tuple[str, ...]) -> str:
        """Helper to format a branch (Observations or Variables) for __repr__."""
        repr_str = f"{indent}├── {branch_name}:\n"
        if cols:
            for i, (col, val) in enumerate(zip(cols, vals, strict=False)):
                connector = "└──" if i == len(cols) - 1 else "├──"
                repr_str += f"{indent}│   {connector} {col}: {val}\n"
        else:
            repr_str += f"{indent}│   └── (all {branch_name.lower()})\n"
        return repr_str

    def __repr__(self):
        """Tree-like representation of the GroupData object."""
        indent = "  "
        repr_str = "GroupByAnnData:\n"

        # Observations branch
        repr_str += self._format_repr_branch(indent, "Observations", self._obs_column_names, self.obs_values)

        # Variables branch
        repr_str += self._format_repr_branch(indent, "Variables", self._var_column_names, self.var_values)

        # AnnData branch
        repr_str += f"{indent}└── AnnData:\n"

        # Get AnnData representation and properly indent it
        adata_repr = self.adata.__repr__()
        adata_lines = adata_repr.split("\n")
        for line in adata_lines:
            repr_str += f"{indent}    {line}\n"

        return repr_str.rstrip()


def _prepare_groups_for_axis(
    adata: ad.AnnData,
    axis_names: pd.Index,
    grouping_func: Callable[[ad.AnnData, Predicates], tuple[Predicates, NwGroupBy]],
    predicates: Predicates | None,
) -> list[tuple[Predicates | None, tuple[str, ...], pd.Index]]:
    """Helper to prepare the list of groups for a given axis (obs or var)."""
    if predicates is not None:
        preds, grouped_result = grouping_func(adata, predicates)
        groups = []
        for group in grouped_result:
            key = tuple(first(group))
            indices = second(group).to_native().index
            groups.append((preds, key, indices))
        return groups
    else:
        return [(None, (), axis_names)]


def _group_by(
    adata: ad.AnnData,
    obs: Predicates,
    var: Predicates,
    copy: bool = False,
) -> ad.AnnData | Generator[GroupByAnndata, None, None]:
    obs_groups = _prepare_groups_for_axis(adata, adata.obs_names, _group_by_obs, obs)
    var_groups = _prepare_groups_for_axis(adata, adata.var_names, _group_by_var, var)

    for (obs_pred, obs_values, obs_idx), (var_pred, var_values, var_idx) in itertools.product(obs_groups, var_groups):
        subset = _construct_adata_from_indices(adata, obs_idx, var_idx)

        # Yield GroupByAnndata with dictionary accessors
        group = GroupByAnndata(
            obs_values=obs_values,
            var_values=var_values,
            adata=subset.copy() if copy else subset,
            _obs_predicates=obs_pred,
            _var_predicates=var_pred,
        )
        yield group
