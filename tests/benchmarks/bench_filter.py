"""Benchmark filter operations on large datasets."""

import annsel as an


def test_filter_obs_100k(large_adata_sparse, benchmark):
    """Benchmark filtering 100K cells by obs column."""
    result = benchmark(lambda: large_adata_sparse.an.filter(obs=an.col("cell_type") == "B cell"))
    assert result.n_obs > 0
    assert result.n_obs < large_adata_sparse.n_obs


def test_filter_obs_quality_100k(large_adata_sparse, benchmark):
    """Benchmark filtering 100K cells by quality threshold."""
    result = benchmark(lambda: large_adata_sparse.an.filter(obs=an.col("quality") > 0.8))
    assert result.n_obs > 0


def test_filter_multicomponent_100k(large_adata_sparse, benchmark):
    """Benchmark multi-component filtering (obs + var)."""
    result = benchmark(
        lambda: large_adata_sparse.an.filter(
            obs=an.col("quality") > 0.7,
            var=an.col("highly_variable") == True,  # noqa: E712
        )
    )
    assert result.n_obs > 0
    assert result.n_vars > 0


def test_select_columns_100k(large_adata_sparse, benchmark):
    """Benchmark column selection on 100K cells."""
    result = benchmark(
        lambda: large_adata_sparse.an.select(obs=an.col("cell_type", "quality"), var=an.col("gene_name"))
    )
    assert len(result.obs.columns) == 2
    assert len(result.var.columns) == 1
