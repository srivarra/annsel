"""Tests for XFrame matrix operations."""

import anndata as ad
import narwhals as nw
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


@pytest.fixture
def adata_with_layers():
    """Create AnnData with X and multiple layers."""
    obs = pd.DataFrame(
        {
            "cell_type": ["B cell", "T cell", "B cell", "T cell", "NK cell"],
            "quality": [0.9, 0.7, 0.8, 0.6, 0.95],
        },
        index=[f"cell_{i}" for i in range(5)],
    )
    var = pd.DataFrame(
        {
            "gene_name": ["GENE1", "GENE2", "GENE3", "GENE4"],
            "highly_variable": [True, False, True, False],
        },
        index=[f"gene_{i}" for i in range(4)],
    )
    # Create X with known values
    X = np.array(
        [
            [10.5, 0.1, 5.2, 1.0],  # cell_0
            [0.2, 2.3, 0.5, 8.1],  # cell_1
            [15.0, 0.3, 8.7, 2.0],  # cell_2
            [0.5, 3.5, 1.1, 12.3],  # cell_3
            [20.0, 0.2, 10.5, 1.5],  # cell_4
        ]
    )

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # Add layers
    adata.layers["raw"] = X * 2  # Raw counts (scaled for testing)
    adata.layers["normalized"] = X / X.sum(axis=1, keepdims=True)  # Normalized

    return adata


@pytest.fixture
def adata_sparse():
    """Create AnnData with sparse X matrix."""
    obs = pd.DataFrame(
        {"cell_type": ["A", "B", "C"]},
        index=[f"cell_{i}" for i in range(3)],
    )
    var = pd.DataFrame(
        {"gene_name": ["G1", "G2", "G3"]},
        index=[f"gene_{i}" for i in range(3)],
    )
    # Create sparse matrix
    X_dense = np.array([[5.0, 0.0, 2.0], [0.0, 3.0, 0.0], [1.0, 0.0, 4.0]])
    X_sparse = sparse.csr_matrix(X_dense)

    return ad.AnnData(X=X_sparse, obs=obs, var=var)


def test_xframe_access(adata_with_layers) -> None:
    """Test accessing X through .nw accessor."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_with_layers)

    # Access X
    x_frame = nwdata.X
    assert x_frame is not None
    assert x_frame.shape == (5, 4)
    assert not x_frame.is_sparse


def test_xframe_layer_access(adata_with_layers) -> None:
    """Test accessing layers through XFrame."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_with_layers)

    # Access raw layer
    raw_layer = nwdata.X.layer("raw")
    assert raw_layer.shape == (5, 4)

    # Access via top-level layer() method
    normalized_layer = nwdata.layer("normalized")
    assert normalized_layer.shape == (5, 4)


def test_xframe_sparse_detection(adata_sparse) -> None:
    """Test sparse matrix detection."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_sparse)

    x_frame = nwdata.X
    assert x_frame.is_sparse


def test_xframe_filter(adata_with_layers) -> None:
    """Test filtering cells by gene expression."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_with_layers)

    # Filter cells where GENE1 > 10
    filtered = nwdata.X.filter(nw.col("gene_0") > 10.0)

    # Should keep cells 0, 2, 4 (GENE1 values: 10.5, 15.0, 20.0)
    assert filtered.n_obs == 3
    assert list(filtered.obs_names) == ["cell_0", "cell_2", "cell_4"]


def test_xframe_select(adata_with_layers) -> None:
    """Test selecting genes via X."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_with_layers)

    # Select first two genes
    selected = nwdata.X.select("gene_0", "gene_1")

    assert selected.n_vars == 2
    assert list(selected.var_names) == ["gene_0", "gene_1"]


def test_xframe_head(adata_with_layers) -> None:
    """Test head operation on X."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_with_layers)

    # Get first 3 cells
    head = nwdata.X.head(3)

    assert head.n_obs == 3
    assert list(head.obs_names) == ["cell_0", "cell_1", "cell_2"]


def test_xframe_sum_mean(adata_sparse) -> None:
    """Test aggregation operations on sparse matrix."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_sparse)

    # Test sum
    total_sum = nwdata.X.sum()
    assert total_sum == 15.0  # 5+2+3+1+4

    # Test sum per gene (axis=0)
    gene_sums = nwdata.X.sum(axis=0)
    np.testing.assert_array_almost_equal(gene_sums, [6.0, 3.0, 6.0])

    # Test sum per cell (axis=1)
    cell_sums = nwdata.X.sum(axis=1)
    np.testing.assert_array_almost_equal(cell_sums, [7.0, 3.0, 5.0])

    # Test mean
    mean_expr = nwdata.X.mean(axis=0)
    np.testing.assert_array_almost_equal(mean_expr, [2.0, 1.0, 2.0])


def test_xframe_layer_filter(adata_with_layers) -> None:
    """Test filtering on specific layer."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_with_layers)

    # Filter on raw layer (X * 2)
    raw_layer = nwdata.layer("raw")
    filtered = raw_layer.filter(nw.col("gene_0") > 20.0)

    # gene_0 in raw: [21, 0.4, 30, 1.0, 40] → cells 0, 2, 4
    assert filtered.n_obs == 3


def test_xframe_to_df_eager(adata_with_layers) -> None:
    """Test DataFrame materialization."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_with_layers)

    # Convert to DataFrame
    df = nwdata.X.to_df()

    assert isinstance(df, pd.DataFrame)
    assert df.shape == (5, 4)
    assert list(df.index) == list(adata_with_layers.obs_names)
    assert list(df.columns) == list(adata_with_layers.var_names)


def test_xframe_layer_error(adata_with_layers) -> None:
    """Test error on non-existent layer."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_with_layers)

    # Try to access non-existent layer
    with pytest.raises(KeyError, match="Layer 'nonexistent' not found"):
        nwdata.X.layer("nonexistent")


def test_xframe_repr(adata_with_layers, adata_sparse) -> None:
    """Test string representation."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)

    # Dense matrix
    nwdata_dense = namespace.from_native(adata_with_layers)
    repr_str = repr(nwdata_dense.X)
    assert "XFrame" in repr_str
    assert "5×4" in repr_str
    assert "dense" in repr_str

    # Sparse matrix
    nwdata_sparse = namespace.from_native(adata_sparse)
    repr_str = repr(nwdata_sparse.X)
    assert "sparse" in repr_str


def test_xframe_integration_with_annsel_api(adata_with_layers) -> None:
    """Test XFrame works with annsel high-level API."""
    # Note: X is not directly on .an (would need to add property)
    # For now, test that with_obs works with X
    combined = adata_with_layers.an.with_obs("cell_type", var_names=["gene_0"])

    assert "cell_type" in combined.to_native().columns
    assert "gene_0" in combined.to_native().columns
    assert len(combined.to_native()) == 5
