"""Utility functions for the AnnData narwhals plugin."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from narwhals.dtypes import DType
    from narwhals.utils import Version


# Placeholder for lit function - in real implementation this would handle literals
def lit(value):
    """Create a literal value."""
    return pd.Series([value])


def evaluate_exprs(df, /, *exprs):
    """Evaluate expressions against a dataframe.

    Parameters
    ----------
    df
        The dataframe context for evaluation.
    exprs
        Expressions to evaluate.

    Returns
    -------
    list[tuple[str, Any]]
        List of (output_name, native_series) tuples.
    """
    native_results: list[tuple[str, Any]] = []
    for expr in exprs:
        native_series_list = expr._call(df)
        output_names = expr._evaluate_output_names(df)
        if expr._alias_output_names is not None:
            output_names = expr._alias_output_names(output_names)
        if len(output_names) != len(native_series_list):  # pragma: no cover
            msg = f"Internal error: got output names {output_names}, but only got {len(native_series_list)} results"
            raise AssertionError(msg)
        native_results.extend(zip(output_names, native_series_list, strict=False))
    return native_results


def native_to_narwhals_dtype(native_dtype: np.dtype, version: Version) -> DType:
    """Convert numpy dtype to narwhals dtype.

    Parameters
    ----------
    native_dtype
        The numpy dtype to convert.
    version
        The narwhals version for compatibility.

    Returns
    -------
    DType
        The corresponding narwhals dtype.
    """
    dtypes = version.dtypes

    # Integer types
    if native_dtype == np.int64:
        return dtypes.Int64()
    if native_dtype == np.int32:
        return dtypes.Int32()
    if native_dtype == np.int16:
        return dtypes.Int16()
    if native_dtype == np.int8:
        return dtypes.Int8()
    if native_dtype == np.uint64:
        return dtypes.UInt64()
    if native_dtype == np.uint32:
        return dtypes.UInt32()
    if native_dtype == np.uint16:
        return dtypes.UInt16()
    if native_dtype == np.uint8:
        return dtypes.UInt8()

    # Float types
    if native_dtype == np.float64:
        return dtypes.Float64()
    if native_dtype == np.float32:
        return dtypes.Float32()

    # Boolean
    if native_dtype == np.bool_:
        return dtypes.Boolean()

    # String/Object
    if native_dtype == np.object_:
        return dtypes.String()

    # Datetime
    if np.issubdtype(native_dtype, np.datetime64):
        return dtypes.Datetime("ns", None)

    # Timedelta
    if np.issubdtype(native_dtype, np.timedelta64):
        return dtypes.Duration("ns")

    return dtypes.Unknown()  # pragma: no cover


def narwhals_to_native_dtype(dtype: DType | type[DType], version: Version) -> np.dtype:
    """Convert narwhals dtype to numpy dtype.

    Parameters
    ----------
    dtype
        The narwhals dtype to convert.
    version
        The narwhals version for compatibility.

    Returns
    -------
    np.dtype
        The corresponding numpy dtype.
    """
    dtypes = version.dtypes

    # Float types
    if dtype == dtypes.Float64:
        return np.dtype("float64")
    if dtype == dtypes.Float32:
        return np.dtype("float32")

    # Integer types
    if dtype == dtypes.Int64:
        return np.dtype("int64")
    if dtype == dtypes.Int32:
        return np.dtype("int32")
    if dtype == dtypes.Int16:
        return np.dtype("int16")
    if dtype == dtypes.Int8:
        return np.dtype("int8")
    if dtype == dtypes.UInt64:
        return np.dtype("uint64")
    if dtype == dtypes.UInt32:
        return np.dtype("uint32")
    if dtype == dtypes.UInt16:
        return np.dtype("uint16")
    if dtype == dtypes.UInt8:
        return np.dtype("uint8")

    # String/Boolean
    if dtype == dtypes.String:
        return np.dtype("object")
    if dtype == dtypes.Boolean:
        return np.dtype("bool")

    # Datetime/Duration
    if dtype == dtypes.Datetime:
        return np.dtype("datetime64[ns]")
    if dtype == dtypes.Duration:
        return np.dtype("timedelta64[ns]")

    msg = f"Unknown dtype: {dtype}"  # pragma: no cover
    raise AssertionError(msg)
