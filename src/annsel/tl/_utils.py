from collections.abc import Callable, Iterable
from functools import reduce
from typing import Any, Self

import anndata as ad
import narwhals as nw
import pandas as pd
from narwhals.typing import DataFrameT, IntoExpr, IntoSeriesT


class PendingOps:
    def __init__(self):
        self._expr: Any = None
        self._pending_ops: list[tuple[Callable, tuple, dict]] = []

    def __getattr__(self, name: str):
        match getattr(self._expr, name, None):
            case None:
                raise AttributeError(f"'{self.__class__.__name__}' has no attribute '{name}'")
            case method:

                def wrapper(*args, **kwargs):
                    self._pending_ops.append((method, args, kwargs))
                    return self

                return wrapper

    # --- binary ---
    def __eq__(self, other: object) -> Self:  # type: ignore[override]
        self._pending_ops.append((self._expr.__eq__, (other,), {}))
        return self

    def __ne__(self, other: object) -> Self:  # type: ignore[override]
        self._pending_ops.append((self._expr.__ne__, (other,), {}))
        return self

    def __and__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__and__, (other,), {}))
        return self

    def __rand__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__rand__, (other,), {}))
        return self

    def __or__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__or__, (other,), {}))
        return self

    def __ror__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__ror__, (other,), {}))
        return self

    def __add__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__add__, (other,), {}))
        return self

    def __radd__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__radd__, (other,), {}))
        return self

    def __sub__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__sub__, (other,), {}))
        return self

    def __rsub__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__rsub__, (other,), {}))
        return self

    def __truediv__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__truediv__, (other,), {}))
        return self

    def __rtruediv__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__rtruediv__, (other,), {}))
        return self

    def __mul__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__mul__, (other,), {}))
        return self

    def __rmul__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__rmul__, (other,), {}))
        return self

    def __le__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__le__, (other,), {}))
        return self

    def __lt__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__lt__, (other,), {}))
        return self

    def __gt__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__gt__, (other,), {}))
        return self

    def __ge__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__ge__, (other,), {}))
        return self

    def __pow__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__pow__, (other,), {}))
        return self

    def __rpow__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__rpow__, (other,), {}))
        return self

    def __floordiv__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__floordiv__, (other,), {}))
        return self

    def __rfloordiv__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__rfloordiv__, (other,), {}))
        return self

    def __mod__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__mod__, (other,), {}))
        return self

    def __rmod__(self, other: Any) -> Self:
        self._pending_ops.append((self._expr.__rmod__, (other,), {}))
        return self

    # --- unary ---
    def __invert__(self) -> Self:
        self._pending_ops.append((self._expr.__invert__, (), {}))
        return self

    def execute(self):
        def apply_op(expr, op_tuple):
            op, args, kwargs = op_tuple
            return op(*args, **kwargs)

        return reduce(apply_op, self._pending_ops, self._expr)


class ObsCol(PendingOps):
    def __init__(self, names: str | Iterable[str]):
        super().__init__()
        self.names = [names] if isinstance(names, str) else list(names)
        self._expr = nw.col(self.names)


class VarCol(PendingOps):
    def __init__(self, names: str | Iterable[str]):
        super().__init__()
        self.names = [names] if isinstance(names, str) else list(names)
        self._expr = nw.col(self.names)


class ObsNames(PendingOps):
    def __init__(self, names: str | Iterable[str]):
        super().__init__()
        self.names = [names] if isinstance(names, str) else list(names)
        self._expr = nw.col(self.names)


class VarNames(PendingOps):
    def __init__(self, names: str | Iterable[str]):
        super().__init__()
        self.names = [names] if isinstance(names, str) else list(names)
        self._expr = nw.col(self.names)


@nw.narwhalify
def _filter_var_names(var_series: IntoSeriesT, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]) -> IntoSeriesT:
    return var_series.filter(*predicates)


@nw.narwhalify
def _filter_obs_names(obs_series: IntoSeriesT, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]) -> IntoSeriesT:
    return obs_series.filter(*predicates)


@nw.narwhalify
def _filter_obs(obs_df: DataFrameT, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]) -> DataFrameT:
    return obs_df.filter(*predicates)


@nw.narwhalify
def _filter_var(var_df: DataFrameT, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]) -> DataFrameT:
    return var_df.filter(*predicates)


def _filter_adata_by_var(adata: ad.AnnData, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]) -> pd.Index:
    return _filter_var(adata.var, *predicates).index


def _filter_adata_by_obs(adata: ad.AnnData, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]) -> pd.Index:
    return _filter_obs(adata.obs, *predicates).index


def _filter_adata_by_var_names(adata: ad.AnnData, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]) -> pd.Index:
    return _filter_var_names(adata.var_names.to_series(), *predicates)


def _filter_adata_by_obs_names(adata: ad.AnnData, *predicates: IntoExpr | Iterable[IntoExpr] | list[bool]) -> pd.Index:
    return _filter_obs_names(adata.obs_names.to_series(), *predicates)
