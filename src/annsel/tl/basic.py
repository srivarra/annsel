from collections.abc import Iterable
from typing import Literal

import anndata as ad
from more_itertools import map_reduce
from narwhals.typing import IntoExpr

from annsel.core.extensions import register_anndata_accessor

from ._utils import (
    ObsCol,
    VarCol,
    _filter_adata_by_obs,
    _filter_adata_by_var,
)


@register_anndata_accessor("an")
class AnnselAccessor:
    """An extension for AnnData."""

    def __init__(self, adata: ad.AnnData):
        self._obj: ad.AnnData = adata

    def filter(self, *predicates: IntoExpr | Iterable[IntoExpr], copy: bool = True) -> ad.AnnData:
        """Filters the AnnData object by the given predicates.

        Parameters
        ----------
        predicates
            The predicates to filter the AnnData object by.

        Returns
        -------
            _description_
        """

        def keyfunc(pred) -> Literal["obs", "var"]:
            match pred:
                case ObsCol():
                    return "obs"
                case VarCol():
                    return "var"
                case _:
                    return "obs"

        def valuefunc(pred):
            match pred:
                case ObsCol() | VarCol():
                    return pred.execute()
                case _:
                    return pred

        grouped_predicates = map_reduce(predicates, keyfunc=keyfunc, valuefunc=valuefunc, reducefunc=list)

        if len(var_predicates := grouped_predicates.get("var", [])) > 0:
            filtered_var_idx = _filter_adata_by_var(self._obj, *var_predicates)
        else:
            filtered_var_idx = self._obj.var_names
        if len(obs_predicates := grouped_predicates.get("obs", [])) > 0:
            filtered_obs_idx = _filter_adata_by_obs(self._obj, *obs_predicates)
        else:
            filtered_obs_idx = self._obj.obs_names
        filtered_adata = self._obj[filtered_obs_idx, filtered_var_idx]
        if copy:
            return filtered_adata.copy()
        else:
            return filtered_adata

    def select(self, *predicates: IntoExpr | Iterable[IntoExpr], copy: bool = True) -> ad.AnnData:
        """#TODO"""
        raise NotImplementedError
        # """Selects the AnnData object by the given predicates."""
        # if len(obs_names_predicates := grouped_predicates.get("obs_names", [])) > 0:
        #     filtered_obs_names = _filter_adata_by_obs_names(self._obj, *obs_names_predicates)
        # else:
        #     filtered_obs_names = []
        # if len(var_names_predicates := grouped_predicates.get("var_names", [])) > 0:
        #     filtered_var_names = _filter_adata_by_var_names(self._obj, *var_names_predicates)
        # else:
        #     filtered_var_names = []
