"""Tests for cross-component operations (obs + X aggregations)."""

import anndata as ad
import narwhals as nw
import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def adata_with_markers():
    """Create AnnData with cell types and marker gene expression."""
    obs = pd.DataFrame(
        {
            "cell_type": ["B cell", "T cell", "B cell", "T cell", "B cell", "T cell"],
            "quality": [0.9, 0.7, 0.85, 0.75, 0.95, 0.65],
            "donor": ["D1", "D1", "D2", "D2", "D1", "D2"],
        },
        index=[f"cell_{i}" for i in range(6)],
    )

    # Marker genes: CD19 (B cell marker), CD3 (T cell marker), CD45 (pan-immune)
    var = pd.DataFrame(
        {
            "gene_name": ["CD19", "CD3", "CD45", "ACTB"],
            "marker_type": ["B cell", "T cell", "Pan", "Housekeeping"],
        },
        index=["CD19", "CD3", "CD45", "ACTB"],
    )

    # Expression matrix with realistic marker patterns
    # B cells: high CD19, low CD3
    # T cells: low CD19, high CD3
    # All: moderate CD45, high ACTB
    X = np.array(
        [
            [15.0, 0.5, 8.0, 20.0],  # B cell - high CD19
            [0.5, 18.0, 7.5, 19.0],  # T cell - high CD3
            [14.0, 0.8, 8.2, 21.0],  # B cell
            [0.8, 17.0, 7.8, 20.5],  # T cell
            [16.0, 0.3, 8.5, 20.5],  # B cell
            [0.3, 16.5, 7.2, 19.5],  # T cell
        ]
    )

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["raw"] = X * 10  # Simulated raw counts

    return adata


def test_with_obs_basic(adata_with_markers) -> None:
    """Test combining obs with expression data."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_with_markers)

    # Combine cell_type from obs with marker gene expression
    combined = nwdata.with_obs(
        "cell_type",
        var_names=["CD19", "CD3"],
    )

    # Verify structure
    df = combined.to_native()
    assert "cell_type" in df.columns
    assert "CD19" in df.columns
    assert "CD3" in df.columns
    assert len(df) == 6


def test_cross_component_aggregation(adata_with_markers) -> None:
    """Test grouping by obs and aggregating X expression."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_with_markers)

    # Combine obs metadata with marker genes
    combined = nwdata.with_obs(
        "cell_type",
        "quality",
        var_names=["CD19", "CD3", "CD45"],
    )

    # Group by cell type and compute statistics
    stats = combined.group_by("cell_type").agg(
        [
            nw.col("quality").mean().alias("mean_quality"),
            nw.col("CD19").mean().alias("CD19_mean"),
            nw.col("CD19").std().alias("CD19_std"),
            nw.col("CD3").mean().alias("CD3_mean"),
            nw.col("CD3").std().alias("CD3_std"),
            nw.col("CD45").mean().alias("CD45_mean"),
            nw.len().alias("n_cells"),
        ]
    )

    result = stats.to_native()

    # Verify results
    assert len(result) == 2  # B cell and T cell
    assert set(result["cell_type"]) == {"B cell", "T cell"}

    # B cells should have high CD19, low CD3
    b_cell_row = result[result["cell_type"] == "B cell"].iloc[0]
    assert b_cell_row["CD19_mean"] > 10  # High CD19
    assert b_cell_row["CD3_mean"] < 2  # Low CD3
    assert b_cell_row["n_cells"] == 3  # 3 B cells

    # T cells should have low CD19, high CD3
    t_cell_row = result[result["cell_type"] == "T cell"].iloc[0]
    assert t_cell_row["CD19_mean"] < 2  # Low CD19
    assert t_cell_row["CD3_mean"] > 15  # High CD3
    assert t_cell_row["n_cells"] == 3  # 3 T cells


def test_with_obs_layer_support(adata_with_markers) -> None:
    """Test with_obs using a specific layer."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_with_markers)

    # Use raw counts layer
    combined = nwdata.with_obs(
        "cell_type",
        var_names=["CD19"],
        layer="raw",
    )

    df = combined.to_native()

    # Raw layer should be X * 10
    assert df["CD19"].iloc[0] > 100  # Should be ~150 (15 * 10)


def test_with_obs_all_obs_cols(adata_with_markers) -> None:
    """Test with_obs without specifying obs_cols (includes all)."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_with_markers)

    # Include all obs columns
    combined = nwdata.with_obs(var_names=["CD19", "CD3"])

    df = combined.to_native()

    # Should have all obs columns
    assert "cell_type" in df.columns
    assert "quality" in df.columns
    assert "donor" in df.columns

    # Plus the selected genes
    assert "CD19" in df.columns
    assert "CD3" in df.columns


def test_marker_expression_by_donor(adata_with_markers) -> None:
    """Test grouping by donor and computing marker expression."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_with_markers)

    # Combine and group by donor
    combined = nwdata.with_obs(
        "donor",
        "cell_type",
        var_names=["CD19", "CD3", "CD45"],
    )

    stats = combined.group_by("donor").agg(
        [
            nw.col("CD19").mean().alias("CD19_mean"),
            nw.col("CD3").mean().alias("CD3_mean"),
            nw.col("CD45").mean().alias("CD45_mean"),
            nw.len().alias("n_cells"),
        ]
    )

    result = stats.to_native()

    assert len(result) == 2  # D1 and D2
    assert all(result["n_cells"] == 3)  # 3 cells per donor


def test_integration_with_annsel_api(adata_with_markers) -> None:
    """Test with_obs through annsel .an accessor."""
    # Use through annsel's .an accessor (unified API)
    combined = adata_with_markers.an.with_obs(
        "cell_type",
        var_names=["CD19", "CD3"],
    )

    # Group and aggregate
    stats = combined.group_by("cell_type").agg(
        [
            nw.col("CD19").mean().alias("CD19_mean"),
            nw.col("CD3").mean().alias("CD3_mean"),
        ]
    )

    result = stats.to_native()

    assert len(result) == 2
    assert "CD19_mean" in result.columns
    assert "CD3_mean" in result.columns
