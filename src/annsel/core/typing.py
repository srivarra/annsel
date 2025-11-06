from collections.abc import Iterable
from typing import TypeAlias

import narwhals as nw

Predicate = nw.Expr | str | nw.Series

# Final recursive type
Predicates: TypeAlias = Predicate | Iterable[Predicate]
