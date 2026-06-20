from collections.abc import Iterable
from typing import Any

import narwhals as nw
from narwhals.group_by import GroupBy, LazyGroupBy

Predicate = nw.Expr | str | nw.Series

# Final recursive type
type Predicates = Predicate | Iterable[Predicate]

# Groupby Types
type NwGroupBy = GroupBy[nw.DataFrame[Any]] | LazyGroupBy[nw.LazyFrame[Any]]
