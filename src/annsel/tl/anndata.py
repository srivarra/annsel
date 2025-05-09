from collections.abc import Callable, Generator, Mapping
from typing import Any, TypeVar

import anndata as ad

from annsel.core.extensions import register_anndata_accessor
from annsel.core.typing import Predicates
from annsel.tl._filter import _filter
from annsel.tl._group_by import GroupByAnndata, _group_by
from annsel.tl._select import _select

# Define type aliases
T = TypeVar("T")


@register_anndata_accessor("an")
class AnnselAccessor:
    """An extension for AnnData."""

    def __init__(self, adata: ad.AnnData):
        self._obj: ad.AnnData = adata

    def filter(
        self,
        obs: Predicates | None = None,
        var: Predicates | None = None,
        x: Predicates | None = None,
        obs_names: Predicates | None = None,
        var_names: Predicates | None = None,
        obsm: Mapping[str, Predicates] | None = None,
        varm: Mapping[str, Predicates] | None = None,
        layer: str | None = None,
        copy: bool = False,
    ) -> ad.AnnData:
        """Filter the AnnData object by the given predicates.

        Parameters
        ----------
        obs
            Predicates to filter the observations by.
        var
            Predicates to filter the variables by.
        x
            Predicates to filter the data by.
        obs_names
            Predicates to filter the observation names by.
        var_names
            Predicates to filter the variable names by.
        obsm
            A mapping of obsm keys and predicates to filter the obsm matrices by.
        varm
            A mapping of varm keys and predicates to filter the varm matrices by.
        layer
            The layer to filter.
        copy
            Whether to return a copy of the AnnData object.

        Returns
        -------
        The filtered AnnData object.

        Examples
        --------
        >>> import annsel as an
        >>> adata = an.datasets.leukemic_bone_marrow_dataset()
        >>> adata.an.filter(
        ...     obs=(
        ...         an.col(["Cell_label"]).is_in(["Classical Monocytes", "CD8+CD103+ tissue resident memory T cells"]),
        ...         an.col(["sex"]) == "male",
        ...     ),
        ...     var=an.col(["vst.mean"]) >= 3,
        ... )
        AnnData object with n_obs × n_vars = 6741 × 458
             obs: 'Cluster_ID', 'donor_id', 'Sample_Tag', 'Cell_label', 'is_primary_data', 'organism_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'assay_ontology_term_id', 'tissue_ontology_term_id', 'Genotype', 'development_stage_ontology_term_id', 'sex_ontology_term_id', 'disease_ontology_term_id', 'cell_type_ontology_term_id', 'suspension_type', 'tissue_type', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid'
             var: 'vst.mean', 'vst.variance', 'vst.variance.expected', 'vst.variance.standardized', 'vst.variable', 'feature_is_filtered', 'Unnamed: 0', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length', 'feature_type'
             uns: 'cell_type_ontology_term_id_colors', 'citation', 'default_embedding', 'schema_reference', 'schema_version', 'title'
             obsm: 'X_bothumap', 'X_pca', 'X_projected', 'X_projectedmean', 'X_tsneni', 'X_umapni'

        """
        _adata = _filter(self._obj, obs, var, x, obs_names, var_names, obsm, varm, layer)
        return _adata if not copy else _adata.copy()

    def select(
        self,
        obs: Predicates | None = None,
        var: Predicates | None = None,
        x: Predicates | None = None,
        copy: bool = False,
    ) -> ad.AnnData:
        """Select the AnnData object by the given predicates.

        Parameters
        ----------
        obs
            Predicates to filter the observations by.
        var
            Predicates to filter the variables by.
        x
            Predicates to filter the data by.
        layer
            The layer to filter.
        copy
            Whether to return a copy of the AnnData object.

        Returns
        -------
        The selected AnnData object.

        Examples
        --------
        >>> import annsel as an
        >>> adata = an.datasets.leukemic_bone_marrow_dataset()
        >>> adata.an.select(
        ...     obs=an.col(["Cell_label", "sex"]),
        ...     var=an.col(["feature_name"]),
        ...     x=an.col(["ENSG00000205336"]),
        ... )
        AnnData object with n_obs × n_vars = 31586 × 1
            obs: 'Cell_label', 'sex'
            var: 'feature_name'
            uns: 'cell_type_ontology_term_id_colors', 'citation', 'default_embedding', 'schema_reference', 'schema_version', 'title'
            obsm: 'X_bothumap', 'X_pca', 'X_projected', 'X_projectedmean', 'X_tsneni', 'X_umapni'
        """
        _adata = _select(self._obj, obs, var, x)
        return _adata if not copy else _adata.copy()

    def group_by(
        self,
        obs: Predicates | None = None,
        var: Predicates | None = None,
        copy: bool = False,
    ) -> ad.AnnData | Generator[GroupByAnndata, None, None]:
        """Group the AnnData object by the given predicates.

        Parameters
        ----------
        obs
            Predicates to group the observations by.
        var
            Predicates to group the variables by.
        copy
            Whether to return a copy of the AnnData object when grouping.
            If obs and var are both None, returns a copy of the original AnnData object if True, otherwise returns the original object.

        Returns
        -------
        The original AnnData object (or a copy) if both `obs` and `var` are `None`.
        Otherwise, returns a generator that yields `GroupByAnndata` objects for each group.

        Examples
        --------
        >>> import annsel as an
        >>> adata = an.datasets.leukemic_bone_marrow_dataset()
        >>> groups = adata.an.group_by(
        ...     obs=an.col(["Cell_label", "sex"]),
        ...     var=an.col(["feature_name"]),
        ... )
        >>> for group in groups:
        ...     print(group.obs_dict, group.var_dict)
        """
        # Handle the case where no grouping is requested
        if obs is None and var is None:
            return self._obj if not copy else self._obj.copy()

        # Otherwise, call the internal generator function
        gb_adata = _group_by(self._obj, obs, var, copy=copy)
        return gb_adata

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

        Examples
        --------
        >>> import annsel as an
        >>> adata = an.datasets.leukemic_bone_marrow_dataset()
        >>> adata.an.pipe(sc.pl.embedding, basis="X_tsneni", color="Cell_label")
        """
        if isinstance(func, tuple):
            func, target = func
            if target in kwargs:
                raise ValueError(f"{target} is both the pipe target and a keyword argument")
            kwargs[target] = self._obj
        else:
            args = (self._obj, *args)
        return func(*args, **kwargs)
