"""Test that annsel API uses narwhals-anndata plugin backend."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import annsel as an


@pytest.fixture
def simple_adata(rng: np.random.Generator) -> ad.AnnData:
    """Create a simple AnnData object for testing."""
    obs = pd.DataFrame(
        {
            "cell_type": ["B cell", "T cell", "B cell", "T cell", "B cell"],
            "quality": [0.9, 0.7, 0.8, 0.6, 0.95],
        },
        index=[f"cell_{i}" for i in range(5)],
    )
    var = pd.DataFrame(
        {
            "gene_name": ["GENE1", "GENE2", "GENE3"],
            "highly_variable": [True, False, True],
        },
        index=[f"gene_{i}" for i in range(3)],
    )
    X = rng.random((5, 3))
    return ad.AnnData(X=X, obs=obs, var=var)


def test_annsel_filter_uses_plugin(simple_adata) -> None:
    """Test that adata.an.filter() uses the narwhals-anndata plugin."""
    # Use the familiar annsel API
    filtered = simple_adata.an.filter(obs=an.col("cell_type") == "B cell")

    # Verify the filtering worked
    assert filtered.n_obs == 3
    assert all(filtered.obs["cell_type"] == "B cell")

    # Verify original data unchanged
    assert simple_adata.n_obs == 5


def test_annsel_select_uses_plugin(simple_adata) -> None:
    """Test that adata.an.select() uses the narwhals-anndata plugin."""
    # Use the familiar annsel API
    selected = simple_adata.an.select(
        obs=an.col("cell_type"),
        var=an.col("gene_name"),
    )

    # Verify selection worked
    assert list(selected.obs.columns) == ["cell_type"]
    assert list(selected.var.columns) == ["gene_name"]

    # Verify original data unchanged
    assert len(simple_adata.obs.columns) == 2
    assert len(simple_adata.var.columns) == 2


def test_with_obs_available(simple_adata) -> None:
    """Test that with_obs() is available on .an accessor."""

    # Test with_obs directly on .an
    combined = simple_adata.an.with_obs("cell_type", var_names=["gene_0"])

    # Verify it returns a LazyFrame
    assert hasattr(combined, "to_native")
    assert "cell_type" in combined.to_native().columns
    assert "gene_0" in combined.to_native().columns


def test_chainable_filter_with_obs(simple_adata) -> None:
    """Test chaining filter with with_obs for aggregation."""
    import narwhals as nw

    # Chain filter → with_obs → agg
    stats = (
        simple_adata.an.filter(obs=an.col("quality") > 0.7)
        .an.with_obs("cell_type", var_names=["gene_0"])
        .group_by("cell_type")
        .agg(
            [
                nw.col("gene_0").mean().alias("gene_0_mean"),
                nw.len().alias("n_cells"),
            ]
        )
    )

    result = stats.to_native()
    assert len(result) >= 1  # At least one cell type
    assert "gene_0_mean" in result.columns
    assert "n_cells" in result.columns
    assert result["n_cells"].sum() == 3  # 3 cells total with quality > 0.7
