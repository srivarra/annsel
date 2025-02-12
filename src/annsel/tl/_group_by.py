import itertools
from collections.abc import Generator, Iterator
from dataclasses import dataclass

import anndata as ad
import narwhals as nw
import pandas as pd
from more_itertools import first, nth
from narwhals.group_by import GroupBy as NwGroupBy
from narwhals.typing import IntoDataFrameT

from annsel.core.typing import GroupBy, Predicates
from annsel.core.utils import _construct_adata_from_indices


@dataclass(frozen=True, slots=True)
class GroupResult:
    """A named tuple that contains the group key and the indices of the group.

    Yields
    ------
        A tuple of the group key and the indices of the group.

    Note
    ----
        Once Narwhals supports `nw.Expr` in `GroupBy` we can remove this, along with the custom AnnselExpr class.
    """

    key: tuple[str, ...]
    indices: pd.Index

    def __iter__(self):
        yield from (self.key, self.indices)


@nw.narwhalify
def _group_by_df(df: IntoDataFrameT, *predicates: Predicates) -> NwGroupBy:
    return df.group_by(*predicates)


def _group_by_obs(adata: ad.AnnData, *predicates: Predicates) -> Generator[GroupResult, None, None]:
    return (
        GroupResult(first(gb), nth(gb, 1).to_native().index)
        for gb in _group_by_df(adata.obs, *(p.names for p in predicates))
    )


def _group_by_var(adata: ad.AnnData, *predicates: Predicates) -> Generator[GroupResult, None, None]:
    return (
        GroupResult(first(gb), nth(gb, 1).to_native().index)
        for gb in _group_by_df(adata.var, *(p.names for p in predicates))
    )


def _group_by(
    adata: ad.AnnData,
    obs: Predicates | None = None,
    var: Predicates | None = None,
    return_group_names: bool = False,
    copy: bool = False,
) -> Iterator[GroupBy]:
    # Create groupings - if no predicate provided, use a single default group.
    obs_gb = list(_group_by_obs(adata, obs)) if obs else [(None, adata.obs_names)]
    var_gb = list(_group_by_var(adata, var)) if var else [(None, adata.var_names)]

    has_obs = obs is not None
    has_var = var is not None

    if not has_obs and not has_var:
        yield adata
        return

    for (obs_groups, obs_idx), (var_groups, var_idx) in itertools.product(obs_gb, var_gb):
        _adata = _construct_adata_from_indices(adata, obs_idx=obs_idx, var_idx=var_idx)
        if copy:
            _adata = _adata.copy()

        if not return_group_names:
            yield _adata
        else:
            match (has_obs, has_var):
                case (False, False):
                    yield (var_groups, _adata)
                case (False, True):
                    yield (obs_groups, _adata)
                case (True, True):
                    yield (obs_groups, var_groups, _adata)
