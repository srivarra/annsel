"""Expression handling for AnnData narwhals plugin."""

from __future__ import annotations

from typing import TYPE_CHECKING

# For now, AnnDataExpr is an alias to narwhals expressions
# since we're delegating expression evaluation to narwhals on obs/var DataFrames
# This will be expanded when we need AnnData-specific expression operations

if TYPE_CHECKING:
    pass

# Type alias for now - expressions are handled by narwhals directly
AnnDataExpr = "nw.Expr"

__all__ = ["AnnDataExpr"]
