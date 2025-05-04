from collections.abc import Iterable
from typing import Any, TypeAlias

import narwhals as nw
from narwhals.group_by import GroupBy, LazyGroupBy

Predicate = nw.Expr | str | nw.Series

# Final recursive type
Predicates: TypeAlias = Predicate | Iterable[Predicate]

# Groupby Types
NwGroupBy: TypeAlias = GroupBy[nw.DataFrame[Any]] | LazyGroupBy[nw.LazyFrame[Any]]
