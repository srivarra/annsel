"""DataFrame wrappers for AnnData narwhals plugin.

This module contains all frame wrappers:
- AnnDataLazyFrame: Top-level wrapper with .obs, .var, .X properties
- ObsFrame: Observation metadata operations
- VarFrame: Variable metadata operations
- XFrame: Expression matrix operations
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal

import anndata as ad
import narwhals as nw
import numpy as np
import pandas as pd
from narwhals._utils import Implementation, ValidateBackendVersion
from narwhals.typing import CompliantLazyFrame

if TYPE_CHECKING:
    from types import ModuleType
    from typing import Self

    from narwhals._utils import _LimitedContext
    from narwhals.dataframe import LazyFrame
    from narwhals.dtypes import DType
    from narwhals.utils import Version
    from scipy import sparse
    from typing_extensions import TypeIs

    from annsel.narwhals_plugin.expr import AnnDataExpr
    from annsel.narwhals_plugin.namespace import AnnDataNamespace


class AnnDataLazyFrame:
    """Main wrapper for AnnData objects exposing obs and var frames.

    This class provides the top-level interface for working with AnnData
    objects in the narwhals ecosystem.
    """

    _implementation = Implementation.UNKNOWN

    def __init__(
        self,
        native_adata: ad.AnnData,
        *,
        version: Version,
        validate_backend_version: bool = False,
    ) -> None:
        """Initialize the AnnData lazy frame.

        Parameters
        ----------
        native_adata
            The native AnnData object to wrap.
        version
            The narwhals version for compatibility.
        validate_backend_version
            Whether to validate backend version compatibility.
        """
        self._native_adata = native_adata
        self._version = version

    @property
    def native(self) -> ad.AnnData:
        """Return the native AnnData object."""
        return self._native_adata

    @property
    def obs(self) -> ObsFrame:
        """Access the observations (obs) frame.

        Returns
        -------
        ObsFrame
            A narwhals-compatible wrapper for the obs DataFrame.
        """
        return ObsFrame(
            self._native_adata,
            axis="obs",
            version=self._version,
        )

    @property
    def var(self) -> VarFrame:
        """Access the variables (var) frame.

        Returns
        -------
        VarFrame
            A narwhals-compatible wrapper for the var DataFrame.
        """
        return VarFrame(
            self._native_adata,
            axis="var",
            version=self._version,
        )

    @property
    def X(self) -> XFrame:
        """Access the expression matrix (X).

        Returns
        -------
        XFrame
            A lazy wrapper for X matrix operations.

        Examples
        --------
        >>> # Filter cells by gene expression
        >>> filtered = nwdata.X.filter(nw.col("GENE1") > 5.0)
        >>>
        >>> # Select specific genes
        >>> selected = nwdata.X.select("GENE1", "GENE2")
        >>>
        >>> # Access layers
        >>> raw_counts = nwdata.X.layer("raw")
        """
        return XFrame(self._native_adata, layer=None, version=self._version)

    def layer(self, name: str) -> XFrame:
        """Access a specific layer.

        Parameters
        ----------
        name
            The layer name to access.

        Returns
        -------
        XFrame
            XFrame wrapper for the specified layer.

        Examples
        --------
        >>> # Access raw counts layer
        >>> raw_layer = nwdata.layer("raw")
        >>> filtered = raw_layer.filter(nw.col("GENE1") > 100)
        """
        return XFrame(self._native_adata, layer=name, version=self._version)

    def to_df(self, layer: str | None = None) -> LazyFrame:
        """Convert X matrix to a DataFrame.

        Parameters
        ----------
        layer
            The layer to convert. If None, uses X.

        Returns
        -------
        LazyFrame
            A narwhals LazyFrame backed by pandas.
        """
        return nw.from_native(self._native_adata.to_df(layer=layer))

    def with_obs(
        self,
        *obs_cols: str,
        var_names: list[str] | None = None,
        layer: str | None = None,
    ) -> LazyFrame:
        """Combine obs metadata with expression values into a single DataFrame.

        This is particularly useful for cross-component aggregations, such as
        grouping by cell type and computing expression statistics for marker genes.

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
            Narwhals LazyFrame with obs columns + expression columns.

        Examples
        --------
        >>> # Combine cell metadata with marker gene expression
        >>> combined = nwdata.with_obs("cell_type", "quality", var_names=["CD3", "CD19", "CD20"], layer="raw")
        >>>
        >>> # Now do cross-component aggregation
        >>> stats = combined.group_by("cell_type").agg(
        ...     [
        ...         nw.col("quality").mean().alias("mean_quality"),
        ...         nw.col("CD3").mean().alias("CD3_mean"),
        ...         nw.col("CD3").std().alias("CD3_std"),
        ...         nw.col("CD19").mean().alias("CD19_mean"),
        ...         nw.col("CD20").mean().alias("CD20_mean"),
        ...     ]
        ... )
        """
        import pandas as pd

        # Get obs columns
        if obs_cols:
            obs_df = self._native_adata.obs[list(obs_cols)]
        else:
            obs_df = self._native_adata.obs

        # Get expression data
        x_df = self._native_adata.to_df(layer=layer)
        if var_names is not None:
            x_df = x_df[var_names]

        # Combine into single DataFrame
        combined = pd.concat([obs_df, x_df], axis=1)

        return nw.from_native(combined)

    def with_var(
        self,
        *var_cols: str,
        obs_names: list[str] | None = None,
        layer: str | None = None,
    ) -> LazyFrame:
        """Combine var metadata with expression values into a single DataFrame.

        This is particularly useful for gene-level analysis, such as grouping by
        gene type and computing expression statistics across cell populations.

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
            Narwhals LazyFrame with var columns + expression columns (transposed).

        Examples
        --------
        >>> # Combine gene metadata with expression across cell types
        >>> combined = nwdata.with_var(
        ...     "gene_type", "highly_variable", obs_names=["B_cell_1", "T_cell_1", "NK_cell_1"], layer="normalized"
        ... )
        >>>
        >>> # Group genes by type and compute expression stats
        >>> stats = combined.group_by("gene_type").agg(
        ...     [
        ...         nw.col("highly_variable").sum().alias("n_variable"),
        ...         nw.col("B_cell_1").mean().alias("mean_in_B"),
        ...         nw.col("T_cell_1").mean().alias("mean_in_T"),
        ...         nw.col("NK_cell_1").mean().alias("mean_in_NK"),
        ...     ]
        ... )
        """
        import pandas as pd

        # Get var columns
        if var_cols:
            var_df = self._native_adata.var[list(var_cols)]
        else:
            var_df = self._native_adata.var

        # Get expression data (transposed: genes as rows, cells as columns)
        x_df = self._native_adata.to_df(layer=layer).T
        if obs_names is not None:
            x_df = x_df[obs_names]

        # Combine into single DataFrame
        combined = pd.concat([var_df, x_df], axis=1)

        return nw.from_native(combined)

    def __narwhals_namespace__(self) -> AnnDataNamespace:
        """Return the narwhals namespace for this backend."""
        from annsel.narwhals_plugin.namespace import AnnDataNamespace

        return AnnDataNamespace(version=self._version)

    def __repr__(self) -> str:
        """Return string representation."""
        return f"AnnDataLazyFrame(n_obs={self._native_adata.n_obs}, n_vars={self._native_adata.n_vars})"


class _BaseAxisFrame(CompliantLazyFrame["AnnDataExpr", "ad.AnnData", "LazyFrame[ad.AnnData]"], ValidateBackendVersion):
    """Abstract base class for ObsFrame and VarFrame.

    This class provides the shared implementation for both obs and var frame operations,
    differing only in axis direction and dataframe access patterns.
    """

    _implementation = Implementation.UNKNOWN

    def __init__(
        self,
        native_adata: ad.AnnData,
        *,
        axis: Literal["obs", "var"],
        version: Version,
        validate_backend_version: bool = False,
    ) -> None:
        """Initialize the axis frame wrapper.

        Parameters
        ----------
        native_adata
            The native AnnData object.
        axis
            Which axis this frame operates on ("obs" or "var").
        version
            The narwhals version for compatibility.
        validate_backend_version
            Whether to validate backend version.
        """
        self._native_adata = native_adata
        self._axis = axis
        self._version = version
        self._cached_schema: dict[str, DType] | None = None
        self._cached_columns: list[str] | None = None
        if validate_backend_version:  # pragma: no cover
            self._validate_backend_version()

    @property
    def _dataframe(self) -> pd.DataFrame:
        """Get the axis-specific DataFrame (obs or var)."""
        return self._native_adata.obs if self._axis == "obs" else self._native_adata.var

    def _slice_adata(self, indices: pd.Index) -> ad.AnnData:
        """Slice AnnData along the appropriate axis."""
        if self._axis == "obs":
            return self._native_adata[indices, :]
        else:
            return self._native_adata[:, indices]

    def _head_slice(self, n: int) -> ad.AnnData:
        """Return first n elements along axis."""
        if self._axis == "obs":
            return self._native_adata[:n, :]
        else:
            return self._native_adata[:, :n]

    @staticmethod
    def _is_native(obj: ad.AnnData | Any) -> TypeIs[ad.AnnData]:
        """Check if object is AnnData."""
        return isinstance(obj, ad.AnnData)

    def to_narwhals(self) -> LazyFrame[ad.AnnData]:
        """Convert to narwhals LazyFrame."""
        return self._version.lazyframe(self, level="lazy")

    def __native_namespace__(self) -> ModuleType:
        """Return the native namespace (anndata module)."""
        import anndata

        return anndata

    def __narwhals_namespace__(self) -> AnnDataNamespace:
        """Return the narwhals namespace."""
        from annsel.narwhals_plugin.namespace import AnnDataNamespace

        return AnnDataNamespace(version=self._version)

    def __narwhals_lazyframe__(self) -> Self:
        """Protocol marker for narwhals lazy frame."""
        return self

    def _with_version(self, version: Version) -> Self:
        """Create new instance with different version."""
        return self.__class__(self._native_adata, axis=self._axis, version=version)

    def _with_native(self, adata: ad.AnnData) -> Self:
        """Create new instance with different AnnData."""
        return self.__class__(adata, axis=self._axis, version=self._version)

    @property
    def native(self) -> ad.AnnData:
        """Return the native AnnData object."""
        return self._native_adata

    @property
    def columns(self) -> list[str]:
        """Return column names."""
        if self._cached_columns is None:
            self._cached_columns = list(self._dataframe.columns)
        return self._cached_columns

    @property
    def schema(self) -> dict[str, DType]:
        """Return schema of DataFrame."""
        if self._cached_schema is None:
            from annsel.narwhals_plugin.utils import native_to_narwhals_dtype

            self._cached_schema = {
                col: native_to_narwhals_dtype(self._dataframe[col].dtype, self._version)
                for col in self._dataframe.columns
            }
        return self._cached_schema

    def filter(self, predicate: AnnDataExpr) -> ad.AnnData:
        """Filter by predicate.

        Parameters
        ----------
        predicate
            The filtering expression.

        Returns
        -------
        ad.AnnData
            A view of the AnnData with filtered elements.
        """
        # Use narwhals on the DataFrame
        nw_df = nw.from_native(self._dataframe)
        filtered_df = nw_df.filter(predicate)

        # Extract the filtered indices
        filtered_indices = filtered_df.to_native().index

        # Return sliced AnnData
        return self._slice_adata(filtered_indices)

    def select(self, *exprs: AnnDataExpr | str) -> ad.AnnData:
        """Select columns.

        Parameters
        ----------
        exprs
            Column names or expressions to select.

        Returns
        -------
        ad.AnnData
            A copy of AnnData with selected columns.
        """
        # Use narwhals on the DataFrame
        nw_df = nw.from_native(self._dataframe)
        selected_df = nw_df.select(*exprs)

        # Create a copy and update the appropriate dataframe
        new_adata = self._native_adata.copy()
        if self._axis == "obs":
            new_adata.obs = selected_df.to_native()
        else:
            new_adata.var = selected_df.to_native()

        return new_adata

    def head(self, n: int = 5) -> ad.AnnData:
        """Return first n elements.

        Parameters
        ----------
        n
            Number of elements to return.

        Returns
        -------
        ad.AnnData
            A view of the first n elements.
        """
        return self._head_slice(n)

    def collect_schema(self) -> dict[str, DType]:
        """Collect schema information."""
        return self.schema

    def with_columns(self, *exprs: AnnDataExpr, **named_exprs: AnnDataExpr) -> ad.AnnData:
        """Add or replace columns in the DataFrame.

        Parameters
        ----------
        exprs
            Expressions to add as new columns.
        named_exprs
            Named expressions where keys become column names.

        Returns
        -------
        ad.AnnData
            Copy of AnnData with modified obs/var columns.

        Examples
        --------
        >>> # Add QC metric to obs
        >>> adata = obs_frame.with_columns(
        ...     pct_mito=(nw.col("mito_counts") / nw.col("total_counts") * 100).alias("pct_mito")
        ... )
        """
        nw_df = nw.from_native(self._dataframe)
        modified_df = nw_df.with_columns(*exprs, **named_exprs)

        new_adata = self._native_adata.copy()
        if self._axis == "obs":
            new_adata.obs = modified_df.to_native()
        else:
            new_adata.var = modified_df.to_native()

        return new_adata

    def drop(self, *columns: str) -> ad.AnnData:
        """Remove columns from the DataFrame.

        Parameters
        ----------
        columns
            Column names to remove.

        Returns
        -------
        ad.AnnData
            Copy of AnnData with specified columns removed.

        Examples
        --------
        >>> # Remove temporary columns
        >>> adata = var_frame.drop("temp_score", "old_annotation")
        """
        nw_df = nw.from_native(self._dataframe)
        modified_df = nw_df.drop(*columns)

        new_adata = self._native_adata.copy()
        if self._axis == "obs":
            new_adata.obs = modified_df.to_native()
        else:
            new_adata.var = modified_df.to_native()

        return new_adata

    def rename(self, mapping: dict[str, str]) -> ad.AnnData:
        """Rename columns in the DataFrame.

        Parameters
        ----------
        mapping
            Dictionary mapping old column names to new names.

        Returns
        -------
        ad.AnnData
            Copy of AnnData with renamed columns.

        Examples
        --------
        >>> # Standardize column names
        >>> adata = var_frame.rename({"gene_name": "feature_name"})
        """
        nw_df = nw.from_native(self._dataframe)
        modified_df = nw_df.rename(mapping)

        new_adata = self._native_adata.copy()
        if self._axis == "obs":
            new_adata.obs = modified_df.to_native()
        else:
            new_adata.var = modified_df.to_native()

        return new_adata

    def sort(
        self,
        by: str | list[str],
        *more_by: str,
        descending: bool | list[bool] = False,
        nulls_last: bool = False,
    ) -> ad.AnnData:
        """Sort rows by column values.

        Parameters
        ----------
        by
            Column name(s) to sort by.
        more_by
            Additional columns to sort by.
        descending
            Sort in descending order.
        nulls_last
            Place null values last.

        Returns
        -------
        ad.AnnData
            View of AnnData with sorted rows.

        Examples
        --------
        >>> # Sort cells by pseudotime
        >>> adata = obs_frame.sort("dpt_pseudotime")
        >>>
        >>> # Sort genes by expression (descending)
        >>> adata = var_frame.sort("mean_expression", descending=True)
        """
        nw_df = nw.from_native(self._dataframe)
        sorted_df = nw_df.sort(by, *more_by, descending=descending, nulls_last=nulls_last)

        # Get sorted indices
        sorted_indices = sorted_df.to_native().index

        return self._slice_adata(sorted_indices)

    def tail(self, n: int = 5) -> ad.AnnData:
        """Return last n elements along axis.

        Parameters
        ----------
        n
            Number of elements to return from end.

        Returns
        -------
        ad.AnnData
            View of last n elements.

        Examples
        --------
        >>> # Get 10 lowest quality cells (after sorting)
        >>> low_quality = obs_frame.sort("n_genes").tail(10)
        """
        if self._axis == "obs":
            return self._native_adata[-n:, :]
        else:
            return self._native_adata[:, -n:]

    def unique(
        self,
        subset: str | list[str] | None = None,
        keep: str = "first",
        maintain_order: bool = False,
    ) -> ad.AnnData:
        """Remove duplicate rows.

        Parameters
        ----------
        subset
            Column(s) to consider for deduplication. If None, all columns.
        keep
            Which duplicate to keep ('first', 'last', 'none', 'any').
        maintain_order
            Preserve original row order.

        Returns
        -------
        ad.AnnData
            View with duplicates removed.

        Examples
        --------
        >>> # Remove duplicate cells by barcode
        >>> adata = obs_frame.unique(subset="cell_barcode")
        """
        nw_df = nw.from_native(self._dataframe)
        unique_df = nw_df.unique(subset=subset, keep=keep, maintain_order=maintain_order)

        # Get unique indices
        unique_indices = unique_df.to_native().index

        return self._slice_adata(unique_indices)


class ObsFrame(_BaseAxisFrame):
    """Narwhals-compatible wrapper for AnnData.obs DataFrame.

    This class enables narwhals operations on the observation metadata,
    maintaining synchronization with the underlying AnnData structure.
    """

    _implementation = Implementation.UNKNOWN

    def __init__(
        self,
        native_adata: ad.AnnData,
        *,
        axis: Literal["obs"],
        version: Version,
        validate_backend_version: bool = False,
    ) -> None:
        """Initialize the obs frame wrapper.

        Parameters
        ----------
        native_adata
            The native AnnData object.
        axis
            Must be "obs" for this class.
        version
            The narwhals version for compatibility.
        validate_backend_version
            Whether to validate backend version.
        """
        super().__init__(
            native_adata,
            axis=axis,
            version=version,
            validate_backend_version=validate_backend_version,
        )

    @classmethod
    def from_native(cls, data: ad.AnnData, /, *, context: _LimitedContext) -> Self:
        """Create from native AnnData object."""
        return cls(data, axis="obs", version=context._version)


class VarFrame(_BaseAxisFrame):
    """Narwhals-compatible wrapper for AnnData.var DataFrame.

    This class enables narwhals operations on the variable metadata,
    maintaining synchronization with the underlying AnnData structure.
    """

    _implementation = Implementation.UNKNOWN

    def __init__(
        self,
        native_adata: ad.AnnData,
        *,
        axis: Literal["var"],
        version: Version,
        validate_backend_version: bool = False,
    ) -> None:
        """Initialize the var frame wrapper.

        Parameters
        ----------
        native_adata
            The native AnnData object.
        axis
            Must be "var" for this class.
        version
            The narwhals version for compatibility.
        validate_backend_version
            Whether to validate backend version.
        """
        super().__init__(
            native_adata,
            axis=axis,
            version=version,
            validate_backend_version=validate_backend_version,
        )

    @classmethod
    def from_native(cls, data: ad.AnnData, /, *, context: _LimitedContext) -> Self:
        """Create from native AnnData object."""
        return cls(data, axis="var", version=context._version)


class XFrame:
    """Lazy wrapper for X matrix operations.

    Provides narwhals-style operations on the expression matrix (X)
    without requiring full materialization to DataFrame.
    """

    def __init__(
        self,
        native_adata: ad.AnnData,
        *,
        layer: str | None = None,
        version: Version,
    ) -> None:
        """Initialize X frame wrapper.

        Parameters
        ----------
        native_adata
            The AnnData object containing the X matrix.
        layer
            The layer to operate on. None for X, string for layers[layer].
        version
            Narwhals version for compatibility.
        """
        self._native_adata = native_adata
        self._layer = layer
        self._version = version

    @property
    def _data(self) -> np.ndarray | sparse.spmatrix:
        """Get the underlying matrix (X or specific layer)."""
        if self._layer is None:
            return self._native_adata.X
        return self._native_adata.layers[self._layer]

    @property
    def shape(self) -> tuple[int, int]:
        """Shape of the matrix (n_obs, n_vars)."""
        return self._data.shape

    @property
    def is_sparse(self) -> bool:
        """Check if the matrix is sparse."""
        from scipy import sparse

        return sparse.issparse(self._data)

    def to_df(self) -> pd.DataFrame:
        """Convert X to DataFrame (eager - loads into memory).

        Returns
        -------
        pd.DataFrame
            DataFrame with obs_names as index, var_names as columns.

        Notes
        -----
        This operation materializes the full matrix into memory,
        which can be expensive for large datasets. Consider filtering
        obs/var first to reduce size.
        """
        return self._native_adata.to_df(layer=self._layer)

    def to_narwhals(self) -> LazyFrame:
        """Convert to narwhals LazyFrame.

        Returns
        -------
        LazyFrame
            Narwhals LazyFrame backed by pandas.
        """
        return nw.from_native(self.to_df())

    def filter(self, predicate: Any) -> ad.AnnData:
        """Filter observations by gene expression values.

        Parameters
        ----------
        predicate
            Narwhals expression to filter rows (observations).
            Expression is evaluated against genes (columns).

        Returns
        -------
        ad.AnnData
            Filtered AnnData view.

        Examples
        --------
        >>> # Filter cells with GENE1 expression > 5
        >>> filtered = nwdata.X.filter(nw.col("GENE1") > 5.0)

        Notes
        -----
        Currently materializes X to DataFrame for filtering.
        Future optimization: Direct sparse matrix filtering.
        """
        # Convert to DataFrame for filtering
        df = self.to_df()
        nw_df = nw.from_native(df)
        filtered = nw_df.filter(predicate)

        # Extract filtered indices
        filtered_indices = filtered.to_native().index

        # Return sliced AnnData
        return self._native_adata[filtered_indices, :]

    def select(self, *exprs: Any) -> ad.AnnData:
        """Select specific genes (variables).

        Parameters
        ----------
        exprs
            Column names or expressions selecting genes.

        Returns
        -------
        ad.AnnData
            AnnData with selected genes.

        Examples
        --------
        >>> # Select specific genes
        >>> selected = nwdata.X.select("GENE1", "GENE2", "GENE3")
        """
        # Convert to DataFrame for selection
        df = self.to_df()
        nw_df = nw.from_native(df)
        selected = nw_df.select(*exprs)

        # Get selected gene names
        selected_genes = selected.to_native().columns

        # Return sliced AnnData
        return self._native_adata[:, selected_genes]

    def head(self, n: int = 5) -> ad.AnnData:
        """Return first n observations.

        Parameters
        ----------
        n
            Number of observations to return.

        Returns
        -------
        ad.AnnData
            View of first n observations.
        """
        return self._native_adata[:n, :]

    def sum(self, axis: int | None = None) -> np.ndarray | float:
        """Sum expression values along axis.

        Parameters
        ----------
        axis
            Axis along which to sum (0=genes, 1=cells, None=all).

        Returns
        -------
        np.ndarray | float
            Sum values.
        """
        if self.is_sparse:
            return np.asarray(self._data.sum(axis=axis)).flatten() if axis is not None else self._data.sum()
        return self._data.sum(axis=axis)

    def mean(self, axis: int | None = None) -> np.ndarray | float:
        """Mean expression values along axis.

        Parameters
        ----------
        axis
            Axis along which to compute mean.

        Returns
        -------
        np.ndarray | float
            Mean values.
        """
        if self.is_sparse:
            return np.asarray(self._data.mean(axis=axis)).flatten() if axis is not None else self._data.mean()
        return self._data.mean(axis=axis)

    def layer(self, name: str) -> XFrame:
        """Access a specific layer.

        Parameters
        ----------
        name
            Layer name to access.

        Returns
        -------
        XFrame
            XFrame wrapper for the specified layer.

        Raises
        ------
        KeyError
            If layer doesn't exist.
        """
        if name not in self._native_adata.layers:
            available = list(self._native_adata.layers.keys())
            msg = f"Layer '{name}' not found. Available layers: {available}"
            raise KeyError(msg)

        return XFrame(self._native_adata, layer=name, version=self._version)

    def __repr__(self) -> str:
        """String representation."""
        layer_str = f"layer='{self._layer}'" if self._layer else "X"
        sparse_str = "sparse" if self.is_sparse else "dense"
        return f"XFrame({layer_str}, {self.shape[0]}×{self.shape[1]}, {sparse_str})"
