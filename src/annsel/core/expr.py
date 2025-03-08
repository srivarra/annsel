from collections.abc import Callable

import narwhals as nw

col: Callable[..., nw.Expr] = nw.col

obs_names: nw.Expr = col("obs_names")
var_names: nw.Expr = col("var_names")
