from collections.abc import Iterator
from typing import Protocol, TypeAlias, runtime_checkable

import anndata as ad
from narwhals.typing import IntoExpr

# Single predicate can be IntoExpr or list[bool]
SinglePredicate: TypeAlias = IntoExpr | list[bool]


@runtime_checkable
class PredicatesCollection(Protocol):
    """Protocol for collections of predicates."""

    def __iter__(self) -> Iterator[SinglePredicate]: ...


# Final recursive type
Predicates: TypeAlias = SinglePredicate | PredicatesCollection


GroupNames: TypeAlias = tuple[str, ...]
GroupBy: TypeAlias = ad.AnnData | tuple[GroupNames, ad.AnnData] | tuple[GroupNames, GroupNames, ad.AnnData]
