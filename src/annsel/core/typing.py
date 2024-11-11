from dataclasses import dataclass
from typing import TypeAlias

import pandas as pd
from narwhals.expr import Expr
from narwhals.series import Series


# from narwhals.typing import IntoExpr
@dataclass
class XIndicies:
    """A dataclass representing the indices of the X matrix of an AnnData object."""

    obs: pd.Index
    var: pd.Index


# Define custom IntoExpr type for documentation mainly as Union["Expr", str, "Series"]
# doesn't play nice with Sphinx autodoc
IntoExpr: TypeAlias = Expr | str | Series
