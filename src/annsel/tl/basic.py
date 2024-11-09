from collections.abc import Callable, Iterable
from functools import reduce
from operator import and_
from typing import Any, TypeVar

import anndata as ad
import pandas as pd
from narwhals.typing import IntoExpr

from annsel.core.extensions import register_anndata_accessor
from annsel.core.utils import _map_predicates
from annsel.tl._filter import (
    FilteredXIndicies,
    _filter_adata_by_obs,
    _filter_adata_by_obs_names,
    _filter_adata_by_var,
    _filter_adata_by_var_names,
    _filter_adata_by_X,
)

T = TypeVar("T")


@register_anndata_accessor("an")
class AnnselAccessor:
    """An extension for AnnData."""

    def __init__(self, adata: ad.AnnData):
        self._obj: ad.AnnData = adata

    def filter(
        self,
        *predicates: IntoExpr | Iterable[IntoExpr],
        copy: bool = True,
        layer: str | None = None,
        keep_sparse: bool = True,
    ) -> ad.AnnData:
        """Filters the AnnData object by the given predicates.

        Parameters
        ----------
        predicates
            The predicates to filter the AnnData object by.

        Returns
        -------
            The filtered AnnData object.

        Examples
        --------
        >>> import annsel as an
        >>> import scanpy as sc
        >>> adata = sc.datasets.pbmc68k_reduced()
        >>> adata.an.filter(
        >>> an.var_names().is_in(["HES4", "TNFRSF4", "SSU72", "PARK7"]),)
        >>>     an.obs_col("bulk_labels").is_in(["CD14+ Monocyte", "CD19+ B"]),
        >>> )
        AnnData object with n_obs × n_vars = 224 × 4
            obs: 'bulk_labels', 'n_genes', 'percent_mito', 'n_counts', 'S_score', 'G2M_score', 'phase', 'louvain'
            var: 'n_counts', 'means', 'dispersions', 'dispersions_norm', 'highly_variable'
            uns: 'bulk_labels_colors', 'louvain', 'louvain_colors', 'neighbors', 'pca', 'rank_genes_groups'
            obsm: 'X_pca', 'X_umap'
            varm: 'PCs'
            obsp: 'distances', 'connectivities'
        """
        grouped_predicates = _map_predicates(*predicates)

        if len(var_predicates := grouped_predicates.var) > 0:
            filtered_var_idx = _filter_adata_by_var(self._obj, *var_predicates)
        else:
            filtered_var_idx = self._obj.var_names

        if len(var_names_predicates := grouped_predicates.var_names) > 0:
            filtered_var_names = _filter_adata_by_var_names(self._obj, *var_names_predicates)
        else:
            filtered_var_names = self._obj.var_names

        if len(obs_predicates := grouped_predicates.obs) > 0:
            filtered_obs_idx = _filter_adata_by_obs(self._obj, *obs_predicates)
        else:
            filtered_obs_idx = self._obj.obs_names

        if len(obs_names_predicates := grouped_predicates.obs_names) > 0:
            filtered_obs_names = _filter_adata_by_obs_names(self._obj, *obs_names_predicates)
        else:
            filtered_obs_names = self._obj.obs_names

        if len(x_predicates := grouped_predicates.x) > 0:
            filtered_x_indices: FilteredXIndicies = _filter_adata_by_X(
                self._obj, *x_predicates, layer=layer, keep_sparse=keep_sparse
            )
        else:
            filtered_x_indices = FilteredXIndicies(self._obj.obs_names, self._obj.var_names)

        final_obs_idx = pd.Series(
            list(reduce(and_, [set(filtered_obs_idx), set(filtered_x_indices.obs), set(filtered_obs_names)]))
        )

        final_var_idx = pd.Series(
            list(reduce(and_, [set(filtered_var_idx), set(filtered_x_indices.var), set(filtered_var_names)]))
        )

        filtered_adata = self._obj[final_obs_idx, final_var_idx]
        if copy:
            return filtered_adata.copy()
        return filtered_adata

    def select(self, *predicates: IntoExpr | Iterable[IntoExpr], copy: bool = True) -> ad.AnnData:
        """Selects the AnnData object by the given predicates.

        Parameters
        ----------
        copy, optional
            _description_, by default True

        Returns
        -------
            _description_
        """
        raise NotImplementedError

    def pipe(self, func: Callable[..., T] | tuple[Callable[..., T], str], *args: Any, **kwargs: Any) -> Any:
        """
        Apply chainable functions ``func(self, *args, **kwargs)`` that expect :class:`~anndata.AnnData`.

        Parameters
        ----------
        func
            Function to apply to this :class:`~anndata.AnnData` object. ``args``, and ``kwargs`` are passed into ``func``.
            Alternatively a ``(callable, data_keyword)`` tuple where ``data_keyword`` is a string indicating
            the keyword of ``callable`` that expects the xarray object.
        *args
            Positional arguments passed into ``func``.
        **kwargs
            Keyword arguments passed into ``func``.

        Returns
        -------
            The return type of ``func``.

        Raises
        ------
        ValueError
            When the data target is both the pipe target and a keyword argument.

        Notes
        -----
        Use ``.pipe`` when chaining together functions that expect AnnData objects, e.g., instead of writing

        .. code:: python

            f(g(h(adata), arg1=a), arg2=b, arg3=cs)

        You can write

        .. code:: python

            (adata.an.pipe(h).pipe(g, arg1=a).an.pipe(f, arg2=b, arg3=c))

        If you have a function that takes the data as (say) the second argument, pass a tuple indicating
        which keyword expects the data. For example suppose ``f`` takes its data as ``arg2``

        .. code:: python

            (adata.an.pipe(h).pipe(g, arg1=a).an.pipe((f, "arg2"), arg1=a, arg3=c))

        See Also
        --------
        pandas.DataFrame.pipe
        xarray.DataArray.pipe
        xarray.Dataset.pipe
        """
        if isinstance(func, tuple):
            func, target = func
            if target in kwargs:
                raise ValueError(f"{target} is both the pipe target and a keyword argument")
            kwargs[target] = self
        else:
            args = (self,) + args
        return func(*args, **kwargs)
