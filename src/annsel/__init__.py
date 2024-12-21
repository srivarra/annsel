from importlib.metadata import version

from . import datasets
from .core import col, obs_names, register_anndata_accessor, var_names
from .tl import AnnselAccessor

__all__ = [
    "AnnselAccessor",
    "col",
    "datasets",
    "obs_names",
    "register_anndata_accessor",
    "var_names",
]

__version__ = version("annsel")
