import itertools
from collections.abc import Generator, Iterator
from dataclasses import dataclass
from typing import Literal

import anndata as ad
import narwhals as nw
import pandas as pd
from more_itertools import first, nth
from narwhals.group_by import GroupBy as NwGroupBy
from narwhals.typing import IntoDataFrameT

from annsel.core.typing import Predicates
from annsel.core.utils import _construct_adata_from_indices

# Define the named tuple at module level


@dataclass(frozen=True, slots=True)
class GroupResult:
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
    sparse_method: Literal["csr", "csc", True, False] | None = None,
) -> Iterator[ad.AnnData | tuple]:
    # Create single-item iterables for None cases
    obs_gb = _group_by_obs(adata, obs) if obs else [(None, adata.obs_names)]
    var_gb = _group_by_var(adata, var) if var else [(None, adata.var_names)]

    # Convert to lists to ensure they can be iterated multiple times
    obs_gb = list(obs_gb)
    var_gb = list(var_gb)

    # If either grouping is empty, yield nothing
    if not obs_gb or not var_gb:
        return

    for (obs_groups, obs_idx), (var_groups, var_idx) in itertools.product(obs_gb, var_gb):
        _adata = _construct_adata_from_indices(adata, obs_idx=obs_idx, var_idx=var_idx, sparse_method=sparse_method)

        if not return_group_names:
            yield _adata
        else:
            match (obs, var):
                case (None, None):
                    yield _adata
                case (None, _):
                    yield (var_groups, _adata)
                case (_, None):
                    yield (obs_groups, _adata)
                case (_, _):
                    yield (obs_groups, var_groups, _adata)
