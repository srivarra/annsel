from collections.abc import Iterator
from typing import Any, Protocol, TypeAlias, runtime_checkable

from narwhals.expr import Expr
from narwhals.series import Series

# Base type for single expressions
IntoExpr: TypeAlias = Expr | str | Series[Any]

# Single predicate can be IntoExpr or list[bool]
SinglePredicate: TypeAlias = IntoExpr | list[bool]


@runtime_checkable
class PredicatesCollection(Protocol):
    """Protocol for collections of predicates."""

    def __iter__(self) -> Iterator[SinglePredicate]: ...


# Final recursive type
Predicates: TypeAlias = SinglePredicate | PredicatesCollection
