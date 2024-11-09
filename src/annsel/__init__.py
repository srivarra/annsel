from importlib.metadata import version

from .core import obs_col, obs_names, register_anndata_accessor, var_col, var_names, x
from .tl import AnnselAccessor

__all__ = ["AnnselAccessor", "obs_col", "var_col", "obs_names", "var_names", "register_anndata_accessor", "x"]

__version__ = version("annsel")
