"""AnnData namespace for narwhals plugin."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from narwhals._compliant.namespace import LazyNamespace
from narwhals._utils import Implementation

if TYPE_CHECKING:
    from collections.abc import Iterable

    import anndata as ad
    from narwhals.dtypes import DType
    from narwhals.typing import ConcatMethod
    from narwhals.utils import Version

    from annsel.narwhals_plugin.dataframe import AnnDataLazyFrame
    from annsel.narwhals_plugin.expr import AnnDataExpr


class AnnDataNamespace(LazyNamespace["AnnDataLazyFrame", "AnnDataExpr", "ad.AnnData"]):
    """Narwhals namespace for AnnData objects.

    This class provides factory methods and top-level functions for working
    with AnnData objects in the narwhals ecosystem.
    """

    _implementation: Implementation = Implementation.UNKNOWN

    def __init__(self, *, version: Version) -> None:
        """Initialize the AnnData namespace.

        Parameters
        ----------
        version
            The narwhals version for compatibility management.
        """
        self._version = version

    def from_native(self, native_object: ad.AnnData) -> AnnDataLazyFrame:
        """Convert a native AnnData object to narwhals-compatible wrapper.

        Parameters
        ----------
        native_object
            The AnnData object to wrap.

        Returns
        -------
        AnnDataLazyFrame
            The wrapped AnnData object.
        """
        from annsel.narwhals_plugin.dataframe import AnnDataLazyFrame

        return AnnDataLazyFrame(native_object, version=self._version)

    @property
    def _expr(self) -> type[AnnDataExpr]:
        """Return the expression class for this backend."""
        from annsel.narwhals_plugin.expr import AnnDataExpr

        return AnnDataExpr

    @property
    def _lazyframe(self) -> type[AnnDataLazyFrame]:
        """Return the lazy frame class for this backend."""
        from annsel.narwhals_plugin.dataframe import AnnDataLazyFrame

        return AnnDataLazyFrame

    def lit(self, value: Any, dtype: DType | type[DType] | None = None) -> AnnDataExpr:
        """Create a literal expression.

        Parameters
        ----------
        value
            The literal value.
        dtype
            Optional dtype for the literal.

        Returns
        -------
        AnnDataExpr
            An expression representing the literal value.
        """
        from annsel.narwhals_plugin.expr import AnnDataExpr
        from annsel.narwhals_plugin.utils import lit, narwhals_to_native_dtype

        def func(_df: AnnDataLazyFrame) -> list:
            if dtype is not None:
                return [lit(value).cast(narwhals_to_native_dtype(dtype, self._version))]
            return [lit(value)]

        return AnnDataExpr(
            func,
            evaluate_output_names=lambda _df: ["literal"],
            alias_output_names=None,
            version=self._version,
        )

    def concat(self, items: Iterable[AnnDataLazyFrame], *, how: ConcatMethod = "vertical") -> AnnDataLazyFrame:
        """Concatenate multiple AnnData objects.

        Parameters
        ----------
        items
            Iterable of AnnData lazy frames to concatenate.
        how
            Concatenation method: "vertical" or "horizontal".

        Returns
        -------
        AnnDataLazyFrame
            Concatenated AnnData object.

        Raises
        ------
        NotImplementedError
            This operation is not yet implemented.
        """
        msg = "concat is not yet implemented for AnnData backend"
        raise NotImplementedError(msg)
