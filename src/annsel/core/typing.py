from collections.abc import Iterator
from typing import Protocol, TypeAlias, runtime_checkable

from narwhals.typing import IntoExpr

# Single predicate can be IntoExpr or list[bool]
SinglePredicate: TypeAlias = IntoExpr | list[bool]


@runtime_checkable
class PredicatesCollection(Protocol):
    """Protocol for collections of predicates."""

    def __iter__(self) -> Iterator[SinglePredicate]: ...


# Final recursive type
Predicates: TypeAlias = SinglePredicate | PredicatesCollection
