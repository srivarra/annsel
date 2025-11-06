"""Narwhals plugin for AnnData objects.

This module implements the narwhals backend protocol for AnnData,
enabling narwhals operations directly on AnnData objects.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData
    from narwhals.utils import Version
    from typing_extensions import TypeIs

    from annsel.narwhals_plugin.namespace import AnnDataNamespace


def __narwhals_namespace__(version: Version) -> AnnDataNamespace:
    """Entry point for narwhals plugin discovery.

    Parameters
    ----------
    version
        The narwhals version for compatibility management.

    Returns
    -------
    AnnDataNamespace
        The namespace object providing access to the plugin's functionality.
    """
    from annsel.narwhals_plugin.namespace import AnnDataNamespace

    return AnnDataNamespace(version=version)


def is_native(native_object: object) -> TypeIs[AnnData]:
    """Type guard to identify AnnData objects.

    Parameters
    ----------
    native_object
        The object to check.

    Returns
    -------
    bool
        True if the object is an AnnData instance.
    """
    try:
        import anndata as ad

        return isinstance(native_object, ad.AnnData)
    except ImportError:
        return False


NATIVE_PACKAGE = "anndata"
"""The name of the native package this plugin wraps."""

__all__ = [
    "NATIVE_PACKAGE",
    "__narwhals_namespace__",
    "is_native",
]
