"""Benchmark cross-component aggregation operations."""

import narwhals as nw
import pytest


@pytest.mark.benchmark
def test_with_obs_small_genes(large_adata_sparse, benchmark):
    """Benchmark with_obs with 10 genes on 100K cells."""
    var_names = large_adata_sparse.var_names[:10].tolist()

    result = benchmark(
        lambda: large_adata_sparse.an.with_obs("cell_type", var_names=var_names)
        .group_by("cell_type")
        .agg([nw.col(var_names[0]).mean(), nw.len()])
        .to_native()
    )
    assert len(result) > 0


@pytest.mark.benchmark
def test_with_obs_many_genes(large_adata_sparse, benchmark):
    """Benchmark with_obs with 100 genes on 100K cells."""
    var_names = large_adata_sparse.var_names[:100].tolist()

    result = benchmark(
        lambda: large_adata_sparse.an.with_obs("cell_type", var_names=var_names)
        .group_by("cell_type")
        .agg([nw.len()])
        .to_native()
    )
    assert len(result) > 0


@pytest.mark.benchmark
def test_with_var_aggregation(large_adata_sparse, benchmark):
    """Benchmark with_var gene-level aggregation."""
    obs_names = large_adata_sparse.obs_names[:1000].tolist()

    result = benchmark(
        lambda: large_adata_sparse.an.with_var("gene_name", obs_names=obs_names)
        .group_by("gene_name")
        .agg([nw.len()])
        .to_native()
    )
    assert len(result) > 0
