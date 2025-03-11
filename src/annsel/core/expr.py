from collections.abc import Iterable

import narwhals as nw


def col(*names: str | Iterable[str]) -> nw.Expr:
    """Create a column expression for the given column name(s) in `obs` or `var`.

    Parameters
    ----------
    column : str or list of str
        The column name(s) to reference

    Returns
    -------
    A column expression that can be used in filter expressions

    Examples
    --------
    >>> import annsel as an
    >>> an.col("gene_1")  # Reference a single column
    >>> an.col(["gene_1", "gene_2"])  # Reference multiple columns
    """
    return nw.col(*names)


obs_names: nw.Expr = col("obs_names")
"""
A column expression representing the observation names `obs_names`.

Examples
--------
>>> import annsel as an
>>> an.obs_names
"""

var_names: nw.Expr = col("var_names")
"""A column expression representing the variable names `var_names`."""
