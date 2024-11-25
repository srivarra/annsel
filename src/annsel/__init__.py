from importlib.metadata import version

from . import datasets
from .core import obs_col, obs_names, register_anndata_accessor, var_col, var_names, x
from .tl import AnnselAccessor

__all__ = [
    "AnnselAccessor",
    "datasets",
    "obs_col",
    "obs_names",
    "register_anndata_accessor",
    "var_col",
    "var_names",
    "x",
]

__version__ = version("annsel")
