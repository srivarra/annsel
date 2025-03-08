import inspect
from collections.abc import Iterable
from functools import reduce
from operator import and_
from typing import TypeVar

import anndata as ad
import narwhals as nw
import pandas as pd
from more_itertools import collapse, nth

from annsel.core.typing import Predicates


def _get_final_indices(
    obj_names: pd.Index,
    *idx: Iterable,
) -> pd.Index:
    # Filter out empty indices
    non_empty_idx = [i for i in idx if len(i) > 0]

    # If all indices were empty, return the original object names
    if not non_empty_idx:
        return obj_names

    # Otherwise, proceed with intersection of non-empty indices
    return obj_names.intersection(pd.Index(list(reduce(and_, map(set, non_empty_idx)))))


def _construct_adata_from_indices(
    adata: ad.AnnData,
    obs_idx: pd.Index,
    var_idx: pd.Index,
) -> ad.AnnData:
    _adata = adata[obs_idx.to_series(), var_idx.to_series()]

    return _adata


T = TypeVar("T")


def second(i: Iterable[T]) -> T:
    """Get the second element from an iterable."""
    snd = nth(i, 1)
    if snd is None:
        raise ValueError("No second element found")
    return snd


def _extract_names_from_expr(*predicates: Predicates) -> tuple[str, ...]:
    names: list[str] = []
    for p in predicates:
        match p:
            case nw.Expr():
                names.extend(collapse(inspect.getclosurevars(p._to_compliant_expr).nonlocals["names"]))
            case str():
                names.append(p)
            case _:
                raise ValueError(f"Unknown predicate type: {type(p)}")
    return tuple(names)
