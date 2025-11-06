from importlib.metadata import version

from . import datasets
from .core import col, obs_names, var_names
from .tl import AnnselAccessor

__all__ = [
    "AnnselAccessor",
    "col",
    "datasets",
    "obs_names",
    "var_names",
]

__version__ = version("annsel")
