import itertools
from collections.abc import Generator, Iterable
from dataclasses import dataclass
from typing import Any

import anndata as ad
import narwhals as nw
from more_itertools import first
from narwhals.typing import Frame

from annsel.core.typing import NwGroupBy, Predicates
from annsel.core.utils import _construct_adata_from_indices, _extract_names_from_expr, second


@nw.narwhalify
def _group_by_df(df: Frame, *names: str) -> NwGroupBy:
    return df.group_by(*names)


def _group_by_obs(adata: ad.AnnData, *predicates: Predicates):
    obs_col_names = _extract_names_from_expr(*predicates)
    _gb = _group_by_df(adata.obs, *obs_col_names)
    return (obs_col_names, _gb)


def _group_by_var(adata: ad.AnnData, *predicates: Predicates):
    var_col_names = _extract_names_from_expr(*predicates)
    _gb = _group_by_df(adata.var, *var_col_names)
    return (var_col_names, _gb)


@dataclass
class GroupByAnndata:
    """Container for grouped AnnData objects and their metadata."""

    # Column names used for grouping
    obs_columns: Iterable[str]
    var_columns: Iterable[str]

    # Values for those columns in this group
    obs_values: tuple[str, ...]
    var_values: tuple[str, ...]

    # The actual data subset
    adata: ad.AnnData

    def obs_dict(self):
        """Dictionary mapping observation column names to their values."""
        return dict(zip(self.obs_columns, self.obs_values, strict=False))

    def var_dict(self):
        """Dictionary mapping variable column names to their values."""
        return dict(zip(self.var_columns, self.var_values, strict=False))

    def __repr__(self):
        """Tree-like representation of the GroupData object."""
        indent = "  "
        repr_str = "GroupByAnnData:\n"

        # Observations branch
        repr_str += f"{indent}├── Observations:\n"
        if self.obs_columns:
            for i, (col, val) in enumerate(zip(self.obs_columns, self.obs_values, strict=False)):
                if i == len(list(self.obs_columns)) - 1:
                    repr_str += f"{indent}│   └── {col}: {val}\n"
                else:
                    repr_str += f"{indent}│   ├── {col}: {val}\n"
        else:
            repr_str += f"{indent}│   └── (all observations)\n"

        # Variables branch
        repr_str += f"{indent}├── Variables:\n"
        if self.var_columns:
            for i, (col, val) in enumerate(zip(self.var_columns, self.var_values, strict=False)):
                if i == len(list(self.var_columns)) - 1:
                    repr_str += f"{indent}│   └── {col}: {val}\n"
                else:
                    repr_str += f"{indent}│   ├── {col}: {val}\n"
        else:
            repr_str += f"{indent}│   └── (all variables)\n"

        # AnnData branch
        repr_str += f"{indent}└── AnnData:\n"

        # Get AnnData representation and properly indent it
        adata_repr = self.adata.__repr__()
        adata_lines = adata_repr.split("\n")
        for line in adata_lines:
            repr_str += f"{indent}    {line}\n"

        return repr_str.rstrip()


def _group_by(
    adata: ad.AnnData,
    obs: Predicates | None = None,
    var: Predicates | None = None,
    copy: bool = False,
) -> Generator[ad.AnnData | GroupByAnndata, Any, None]:
    if obs is None and var is None:
        yield adata
        return

    if obs is not None:
        obs_col_names, obs_grouped = _group_by_obs(adata, obs)
        obs_groups = []

        # Extract groups from the grouped object
        for group in obs_grouped:
            key = tuple(first(group))  # Get group values
            indices = second(group).to_native().index  # Get indices
            obs_groups.append((obs_col_names, key, indices))
    else:
        # Default: one group with all observations
        obs_groups = [([], (), adata.obs_names)]

    if var is not None:
        var_col_names, var_grouped = _group_by_var(adata, var)
        var_groups = []

        # Extract groups from the grouped object
        for group in var_grouped:
            key = tuple(first(group))  # Get group values
            indices = second(group).to_native().index  # Get indices
            var_groups.append((var_col_names, key, indices))
    else:
        # Default: one group with all variables
        var_groups = [([], (), adata.var_names)]

    for (obs_cols, obs_values, obs_idx), (var_cols, var_values, var_idx) in itertools.product(obs_groups, var_groups):
        subset = _construct_adata_from_indices(adata, obs_idx, var_idx)
        if copy:
            subset = subset.copy()

        # Yield GroupData with dictionary accessors
        yield GroupByAnndata(
            obs_columns=obs_cols, var_columns=var_cols, obs_values=obs_values, var_values=var_values, adata=subset
        )
