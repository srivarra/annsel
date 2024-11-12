from abc import ABC, abstractmethod
from collections.abc import Callable, Iterable
from functools import wraps
from typing import Concatenate, ParamSpec, TypeVar

import anndata as ad
import pandas as pd
from narwhals.typing import IntoExpr

from annsel.core.typing import XIndicies
from annsel.core.utils import GroupedPredicates, _map_predicates

T = TypeVar("T")
P = ParamSpec("P")


def predicate_guard(
    predicate_type: str,
) -> Callable[[Callable[P, T]], Callable[Concatenate["NarwhalsMethod", P], T | None]]:
    """Decorator that only executes the method if corresponding predicates exist."""

    def decorator(func: Callable[..., T]) -> Callable[..., T | None]:
        @wraps(func)
        def wrapper(self: NarwhalsMethod, *args, **kwargs) -> T | None:
            predicates = getattr(self._predicates, predicate_type)
            if len(predicates) > 0:
                return func(self, *args, **kwargs)
            return None

        return wrapper

    return decorator


class NarwhalsMethod(ABC):
    """Base class for Narwhals methods that operate on AnnData objects.

    This abstract base class defines the interface for methods that filter or transform
    AnnData objects using predicates. Subclasses must implement the abstract methods
    to define specific filtering/transformation behavior.

    Parameters
    ----------
    adata
        The AnnData object to operate on
    predicates
        The predicates to apply to the AnnData object
    """

    def __init__(self, adata: ad.AnnData, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]):
        self._adata = adata
        self._predicates: GroupedPredicates = _map_predicates(*predicates)

    @abstractmethod
    def __call__(self, *args, **kwargs):
        """Execute the command on the AnnData object."""
        raise NotImplementedError

    @abstractmethod
    def _run_var_predicates(self) -> pd.Index | None:
        raise NotImplementedError

    @abstractmethod
    def _run_var_names_predicates(self) -> pd.Index | None:
        raise NotImplementedError

    @abstractmethod
    def _run_obs_predicates(self) -> pd.Index | None:
        raise NotImplementedError

    @abstractmethod
    def _run_obs_names_predicates(self) -> pd.Index | None:
        raise NotImplementedError

    @abstractmethod
    def _run_x_predicates(self, layer: str | None = None, keep_sparse: bool = True) -> XIndicies | None:
        raise NotImplementedError

    @abstractmethod
    def _finalize_indices_obs(self, *idx: pd.Index) -> pd.Index:
        raise NotImplementedError

    @abstractmethod
    def _finalize_indices_var(self, *idx: pd.Index) -> pd.Index:
        raise NotImplementedError
