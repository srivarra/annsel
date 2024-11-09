from collections.abc import Callable, Iterable
from typing import Any

from narwhals.utils import flatten

from .expr import ObsExpr, ObsNamesExpr, VarExpr, VarNamesExpr, XExpr


def _with_names(func: Callable) -> Callable:
    def wrapper(plx: Any, *names: str | Iterable[str]) -> Any:
        return plx.col(*flatten(names))

    return wrapper


@_with_names
def _func(plx: Any, *names: str | Iterable[str]) -> Any:
    return plx.col(*names)


class ObsCol:
    """Select columns from the :obj:`~anndata.AnnData.obs` DataFrame of an :obj:`~anndata.AnnData` object."""

    def __call__(self, *names: str | Iterable[str]) -> ObsExpr:
        """Select columns from the :obj:`~anndata.AnnData.obs` DataFrame of an :obj:`~anndata.AnnData` object.

        This is a wrapper around the `narwhals.col` function


        Parameters
        ----------
        names
            The names of the obs columns to select.

        Returns
        -------
        A `narwhals.Expr` object representing the selected columns.
        """
        return ObsExpr(lambda plx: _func(plx, *names))


class VarCol:
    """Select columns from the :obj:`~anndata.AnnData.var` DataFrame of an :obj:`~anndata.AnnData` object."""

    def __call__(self, *names: str | Iterable[str]) -> VarExpr:
        """Select columns from the :obj:`~anndata.AnnData.var` DataFrame of an :obj:`~anndata.AnnData` object.

        This is a wrapper around the `narwhals.col` function

        Parameters
        ----------
        names
            The names of the var columns to select.

        Returns
        -------
        An `annsel.core.expr.VarExpr` object representing the selected columns.
        """
        return VarExpr(lambda plx: _func(plx, *names))


class ObsNames:
    """Selects the :obj:`~anndata.AnnData.obs_names` index of the :obj:`~anndata.AnnData` object."""

    def __call__(self) -> ObsExpr:
        """Selects the :obj:`~anndata.AnnData.obs_names` index of the :obj:`~anndata.AnnData` object.

        Returns
        -------
        An `annsel.core.expr.ObsExpr` object representing the selected columns.
        """
        return ObsNamesExpr(lambda plx: _func(plx, "obs_names"))


class VarNames:
    """Selects the :obj:`~anndata.AnnData.var_names` index of the :obj:`~anndata.AnnData` object."""

    def __call__(self) -> VarExpr:
        """Selects the :obj:`~anndata.AnnData.var_names` index of the :obj:`~anndata.AnnData` object.

        Returns
        -------
        An `annsel.core.expr.VarExpr` object representing the selected columns.
        """
        return VarNamesExpr(lambda plx: _func(plx, "var_names"))


class X:
    """Selects the :obj:`~anndata.AnnData.X` matrix of the :obj:`~anndata.AnnData` object."""

    def __call__(self, names: str | Iterable[str]) -> XExpr:
        """Selects the :obj:`~anndata.AnnData.X` matrix of the :obj:`~anndata.AnnData` object.

        Returns
        -------
        An `annsel.core.expr.XExpr` object representing the selected columns.
        """
        return XExpr(lambda plx: _func(plx, *names))


obs_col: ObsCol = ObsCol()
var_col: VarCol = VarCol()
obs_names: ObsNames = ObsNames()
var_names: VarNames = VarNames()
x: X = X()
