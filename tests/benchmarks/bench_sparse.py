"""Benchmark sparse matrix operations."""

import pytest

import annsel as an


@pytest.mark.benchmark
def test_sparse_filter_preservation(large_adata_sparse, benchmark):
    """Benchmark filtering sparse 100K×1K matrix."""

    def filter_sparse():
        result = large_adata_sparse.an.filter(obs=an.col("n_genes") > 1000)
        # Verify sparsity preserved
        assert hasattr(result.X, "nnz")  # Still sparse
        return result

    result = benchmark(filter_sparse)
    assert result.n_obs > 0


@pytest.mark.benchmark
def test_sparse_multicomponent_filter(large_adata_sparse, benchmark):
    """Benchmark multi-component filtering on sparse data."""
    result = benchmark(
        lambda: large_adata_sparse.an.filter(obs=an.col("quality") > 0.5, var=an.col("mean_expression") > 1.0)
    )
    assert result.n_obs > 0
    assert result.n_vars > 0
