from collections.abc import Callable
from typing import Any, TypeVar

import anndata as ad
from anndata import register_anndata_namespace

from annsel.core.typing import Predicates
from annsel.tl._filter import _filter_plugin
from annsel.tl._select import _select_plugin

# Define type aliases
T = TypeVar("T")


@register_anndata_namespace("an")
class AnnselAccessor:
    """An extension for AnnData.

    This accessor provides a convenient API for DataFrame-style operations
    on AnnData objects, powered by the narwhals-anndata plugin.
    """

    def __init__(self, adata: ad.AnnData) -> None:
        self._obj: ad.AnnData = adata

    def filter(
        self,
        obs: Predicates | None = None,
        var: Predicates | None = None,
        x: Predicates | None = None,
        obs_names: Predicates | None = None,
        var_names: Predicates | None = None,
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
        layer
            The layer to filter (for x filtering).
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
        # Use the new plugin-based filter implementation
        _adata = _filter_plugin(self._obj, obs, var, x, obs_names, var_names, layer)
        return _adata if not copy else _adata.copy()

    def select(
        self,
        obs: Predicates | None = None,
        var: Predicates | None = None,
        x: Predicates | None = None,
        layer: str | None = None,
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
        # Use the new plugin-based select implementation
        _adata = _select_plugin(self._obj, obs, var, x, layer)
        return _adata if not copy else _adata.copy()

    def with_obs(
        self,
        *obs_cols: str,
        var_names: list[str] | None = None,
        layer: str | None = None,
    ):
        """Combine obs metadata with expression values for cross-component aggregation.

        This method enables powerful workflows like grouping by cell type and
        computing marker gene expression statistics.

        Parameters
        ----------
        obs_cols
            Obs column names to include. If empty, includes all obs columns.
        var_names
            Gene names to include from X. If None, includes all genes.
        layer
            Layer to use for expression values. If None, uses X.

        Returns
        -------
        LazyFrame
            Narwhals LazyFrame ready for group_by().agg() operations.

        Examples
        --------
        >>> import annsel as an
        >>> import narwhals as nw
        >>> # Group by cell type, compute marker expression stats
        >>> stats = (
        ...     adata.an.with_obs("cell_type", var_names=["CD3", "CD19", "CD20"], layer="raw")
        ...     .group_by("cell_type")
        ...     .agg(
        ...         [
        ...             nw.col("CD3").mean().alias("CD3_mean"),
        ...             nw.col("CD19").mean().alias("CD19_mean"),
        ...         ]
        ...     )
        ... )
        """
        from narwhals._utils import Version

        from annsel.narwhals_plugin import __narwhals_namespace__

        namespace = __narwhals_namespace__(Version.MAIN)
        lazy_frame = namespace.from_native(self._obj)
        return lazy_frame.with_obs(*obs_cols, var_names=var_names, layer=layer)

    def with_var(
        self,
        *var_cols: str,
        obs_names: list[str] | None = None,
        layer: str | None = None,
    ):
        """Combine var metadata with expression values for gene-level aggregation.

        Parameters
        ----------
        var_cols
            Var column names to include. If empty, includes all var columns.
        obs_names
            Cell names to include from X. If None, includes all cells.
        layer
            Layer to use for expression values. If None, uses X.

        Returns
        -------
        LazyFrame
            Narwhals LazyFrame ready for group_by().agg() operations.

        Examples
        --------
        >>> import annsel as an
        >>> import narwhals as nw
        >>> # Group genes by type, compute expression stats
        >>> stats = (
        ...     adata.an.with_var("gene_type", obs_names=["B_cell_1", "T_cell_1"])
        ...     .group_by("gene_type")
        ...     .agg(
        ...         [
        ...             nw.col("B_cell_1").mean(),
        ...             nw.col("T_cell_1").mean(),
        ...         ]
        ...     )
        ... )
        """
        from narwhals._utils import Version

        from annsel.narwhals_plugin import __narwhals_namespace__

        namespace = __narwhals_namespace__(Version.MAIN)
        lazy_frame = namespace.from_native(self._obj)
        return lazy_frame.with_var(*var_cols, obs_names=obs_names, layer=layer)

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
                msg = f"{target} is both the pipe target and a keyword argument"
                raise ValueError(msg)
            kwargs[target] = self._obj
        else:
            args = (self._obj, *args)
        return func(*args, **kwargs)
