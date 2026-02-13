from collections.abc import Iterable
from functools import reduce
from operator import and_
from typing import TypeVar

import anndata as ad
import narwhals as nw
import pandas as pd
from more_itertools import nth
from narwhals._expression_parsing import ExprKind

from annsel.core.typing import Predicates

T = TypeVar("T")


def _get_final_indices(
    obj_names: pd.Index,
    *idx: Iterable,
) -> pd.Index:
    """Given an iterable of indices, return the intersection of the indices w.r.t. the var_names or the obs_names.

    If all indices are empty, return an empty index.

    Parameters
    ----------
    obj_names
        The names of the objects to filter.
    idx
        The indices to filter.

    Returns
    -------
    The intersection of the indices.
    """
    # Filter out empty indices
    non_empty_idx = [i for i in idx if len(i) > 0]

    # If all indices were empty, return an empty index
    if not non_empty_idx:
        return pd.Index([])

    # Otherwise, compute the intersection of non-empty indices
    return obj_names.intersection(pd.Index(list(reduce(and_, map(set, non_empty_idx)))))


def _construct_adata_from_indices(
    adata: ad.AnnData,
    obs_idx: pd.Index,
    var_idx: pd.Index,
) -> ad.AnnData:
    _adata = adata[obs_idx.to_series(), var_idx.to_series()]

    return _adata


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
                # doesn't handle nested expressions. E.g. for
                # (nw.col('x') > 5) & (nw.col('y') < 10) it gives ('x',)
                for node in p._nodes:
                    if node.kind == ExprKind.COL:
                        names.extend(node.kwargs["names"])
            case str():
                names.append(p)
            case _:
                raise ValueError(f"Unknown predicate type: {type(p)}")
    return tuple(names)
