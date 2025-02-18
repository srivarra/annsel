from collections.abc import Callable, Iterable
from typing import Any, Protocol

import narwhals as nw
from narwhals._expression_parsing import ExprKind, ExprMetadata
from narwhals.utils import flatten


class AnnselExpr(nw.Expr):
    """A wrapper for the `narwhals.Expr` class."""

    def __init__(
        self,
        to_compliant_expr: Callable[[Any], Any],
        metadata: ExprMetadata,
        *names: Iterable[str],
    ) -> None:
        self._to_compliant_expr = to_compliant_expr
        self._metadata = metadata

        self.names = names


class ColProtocol(Protocol):
    """Protocol for column selection."""

    def __call__(self, *names: str | Iterable[str]) -> AnnselExpr:
        """Protocol for column selection.

        Parameters
        ----------
        *names
            Names of columns to select. Can be strings or iterables of strings.

        Returns
        -------
        AnnselExpr
            A wrapper around the narwhals.Expr class for column selection.
        """
        ...


class Col:
    """A wrapper around the `narwhals.col` function."""

    names: Iterable[str]

    def __call__(self, *names: str | Iterable[str]) -> AnnselExpr:
        """
        Select columns from a DataFrame component of an AnnData object.

        Parameters
        ----------
        *names
            Names of columns to select. Can be strings or iterables of strings.

        Returns
        -------
        A wrapper around the `narwhals.Expr` class.
        """

        def column_selector(plx: Any):
            return plx.col(*flatten(names))

        return AnnselExpr(
            column_selector, ExprMetadata(kind=ExprKind.TRANSFORM, is_order_dependent=False), *flatten(names)
        )


col: Col = Col()
obs_names: AnnselExpr = col("obs_names")
var_names: AnnselExpr = col("var_names")
