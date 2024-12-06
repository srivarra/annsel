from collections.abc import Callable, Iterable
from functools import wraps
from inspect import signature
from typing import Any

from narwhals.utils import flatten

from annsel.core.expr import ObsExpr, ObsNamesExpr, VarExpr, VarNamesExpr, XExpr


def _with_names(func: Callable) -> Callable:
    sig = signature(func)

    @wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        # Bind the arguments to the function signature
        bound_args = sig.bind(*args, **kwargs)
        bound_args.apply_defaults()

        # Extract plx and names from bound arguments
        plx = bound_args.arguments["plx"]
        names = bound_args.arguments["names"]

        # Flatten the names and call plx.col
        return plx.col(*flatten(names))

    return wrapper


@_with_names
def _func(plx: Any, *names: str | Iterable[str]) -> Any:
    return plx.col(*names)


class ObsCol:
    """Select columns from the :obj:`~anndata.AnnData.obs` DataFrame of an :obj:`~anndata.AnnData` object.

    An instance of this class is exported under the name `obs_col`. It can be used as
    though it were a function by calling, for example, `an.obs_col("foo")`.
    See the :func:`__call__` method for further documentation.

    """

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

        Examples
        --------
        >>> import annsel as an
        >>> an.obs_col(["cell_cluster"]).is_in(["Tumor", "Stroma"])  # Select all obs_col > 100
        """
        return ObsExpr(lambda plx: _func(plx, *names))


class VarCol:
    """Select columns from the :obj:`~anndata.AnnData.var` DataFrame of an :obj:`~anndata.AnnData` object.

    An instance of this class is exported under the name `obs_col`. It can be used as
    though it were a function by calling, for example, `an.obs_col("foo")`.
    See the :func:`__call__` method for further documentation.
    """

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

        Examples
        --------
        >>> import annsel as an
        >>> an.var_col(["mean"]) > 1  # Select all var_col with "mean" > 1
        """
        return VarExpr(lambda plx: _func(plx, *names))


class ObsNames:
    """Selects the :obj:`~anndata.AnnData.obs_names` index of the :obj:`~anndata.AnnData` object."""

    def __call__(self) -> ObsNamesExpr:
        """Selects the :obj:`~anndata.AnnData.obs_names` index of the :obj:`~anndata.AnnData` object.

        Returns
        -------
        An `annsel.core.expr.ObsExpr` object representing the selected columns.

        Examples
        --------
        >>> import annsel as an
        >>> an.obs_names(["obs1", "obs2", "obs3"])  # select obs_names
        >>> an.obs_names().ends_with("100")  # select all obs_names ending with "100"
        """
        return ObsNamesExpr(lambda plx: _func(plx, "obs_names"))


class VarNames:
    """Selects the :obj:`~anndata.AnnData.var_names` index of the :obj:`~anndata.AnnData` object."""

    def __call__(self) -> VarNamesExpr:
        """Selects the :obj:`~anndata.AnnData.var_names` index of the :obj:`~anndata.AnnData` object.

        Returns
        -------
        An `annsel.core.expr.VarExpr` object representing the selected columns.

        Examples
        --------
        >>> import annsel as an
        >>> an.var_names(["var1", "var2", "var3"])  # select var_names
        >>> an.var_names().ends_with("100")  # select all var_names ending with "100"
        """
        return VarNamesExpr(lambda plx: _func(plx, "var_names"))


class X:
    """Selects the :obj:`~anndata.AnnData.X` matrix of the :obj:`~anndata.AnnData` object."""

    def __call__(self, names: str | Iterable[str]) -> XExpr:
        """Selects the :obj:`~anndata.AnnData.X` matrix of the :obj:`~anndata.AnnData` object.

        Returns
        -------
        An `annsel.core.expr.XExpr` object representing the selected columns.

        Examples
        --------
        >>> import annsel as an
        >>> an.x(["gene1", "gene2", "gene3"]) > 100  # Select all x columns with expression > 100
        """
        return XExpr(lambda plx: _func(plx, *names))


obs_col: ObsCol = ObsCol()
var_col: VarCol = VarCol()
obs_names: ObsNames = ObsNames()
var_names: VarNames = VarNames()
x: X = X()
