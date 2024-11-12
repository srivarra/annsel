from collections.abc import Callable, Iterable
from typing import Any, Literal, TypeVar

import anndata as ad

from annsel.core.extensions import register_anndata_accessor
from annsel.core.typing import IntoExpr
from annsel.core.utils import _handle_sparse_method
from annsel.tl._filter import (
    FilterAnnData,
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
        layer: str | None = None,
        keep_sparse: bool = True,
        sparse_method: Literal["csr", "csc"] | None = None,
        copy: bool = False,
    ) -> ad.AnnData:
        """Filters the AnnData object by the given predicates.

        Parameters
        ----------
        predicates
            The predicates to filter the AnnData object by.
        layer
            The layer to filter the AnnData object by.
        keep_sparse
            Whether to keep the sparse matrix.
        sparse_method
            Convert X to a sparse array if desired. `"csr"` and `"csc"` are supported.
        copy
            Whether to return a copy of the AnnData object, or return a view of the original.
            Defaults to `False`, which returns a view.

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
        filter_adata = FilterAnnData(self._obj, *predicates)
        final_obs_idx, final_var_idx = filter_adata(layer=layer, keep_sparse=keep_sparse, sparse_method=sparse_method)

        filtered_adata = self._obj[final_obs_idx, final_var_idx]

        filtered_adata = _handle_sparse_method(filtered_adata, sparse_method)
        if copy:
            return filtered_adata.copy()
        return filtered_adata

    def select(self, *predicates: IntoExpr | Iterable[IntoExpr], copy: bool = True) -> ad.AnnData:
        """Selects the AnnData object by the given predicates.

        Parameters
        ----------
        copy
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
            kwargs[target] = self._obj
        else:
            args = (self._obj, *args)
        return func(*args, **kwargs)
