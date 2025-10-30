from collections.abc import Iterable
from typing import TypeVar

import anndata as ad
import narwhals as nw
import pandas as pd
from more_itertools import nth

from annsel.core.typing import Predicates

T = TypeVar("T")


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
        msg = "No second element found"
        raise ValueError(msg)
    return snd


def _extract_names_from_expr(*predicates: Predicates) -> tuple[str, ...]:
    """Extract column names from narwhals expressions or strings.

    Works with current narwhals version by accessing the _nodes structure.

    Parameters
    ----------
    predicates
        Narwhals expressions or column name strings.

    Returns
    -------
    tuple[str, ...]
        Extracted column names.
    """
    names: list[str] = []
    for p in predicates:
        match p:
            case nw.Expr():
                # Extract from narwhals expression nodes
                for node in p._nodes:
                    if node.kind.name == "COL" and "names" in node.kwargs:
                        names.extend(node.kwargs["names"])
            case str():
                names.append(p)
            case _:
                msg = f"Unknown predicate type: {type(p)}"
                raise ValueError(msg)
    return tuple(names)
