from collections.abc import Callable
from typing import Any

import narwhals as nw


class ObsExpr(nw.Expr):
    """A Obs DataFrame wrapper for the `narwhals.Expr` class."""

    def __init__(self, call: Callable[[Any], Any]) -> None:
        super().__init__(call)


class VarExpr(nw.Expr):
    """A Var DataFrame wrapper for the `narwhals.Expr` class."""

    def __init__(self, call: Callable[[Any], Any]) -> None:
        super().__init__(call)


class XExpr(nw.Expr):
    """A X DataFrame wrapper for the `narwhals.Expr` class."""

    def __init__(self, call: Callable[[Any], Any]) -> None:
        super().__init__(call)


class ObsNamesExpr(nw.Expr):
    """A Obs Names DataFrame wrapper for the `narwhals.Expr` class."""

    def __init__(self, call: Callable[[Any], Any]) -> None:
        super().__init__(call)


class VarNamesExpr(nw.Expr):
    """A Var Names DataFrame wrapper for the `narwhals.Expr` class."""

    def __init__(self, call: Callable[[Any], Any]) -> None:
        super().__init__(call)
