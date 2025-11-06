"""Fixtures for performance benchmarks."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


@pytest.fixture(scope="session")
def large_adata_sparse():
    """Create large sparse AnnData: 100K cells × 1K genes.

    Realistic single-cell RNA-seq scale with:
    - Cell type annotations
    - QC metrics
    - Sparse expression matrix (10% density)
    - One layer for testing
    """
    # Use fixed seed for reproducible benchmarks
    rng = np.random.default_rng(42)

    n_obs = 100_000
    n_vars = 1_000

    # Obs metadata
    obs = pd.DataFrame(
        {
            "cell_type": rng.choice(["B cell", "T cell", "NK cell", "Monocyte"], n_obs),
            "quality": rng.random(n_obs),
            "n_genes": rng.integers(200, 5000, n_obs),
            "donor": rng.choice(["D1", "D2", "D3"], n_obs),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )

    # Var metadata
    var = pd.DataFrame(
        {
            "gene_name": [f"GENE{i}" for i in range(n_vars)],
            "highly_variable": rng.random(n_vars) > 0.8,
            "mean_expression": rng.random(n_vars) * 10,
        },
        index=[f"gene_{i}" for i in range(n_vars)],
    )

    # Sparse expression matrix (10% density - realistic for scRNA-seq)
    X = sparse.random(n_obs, n_vars, density=0.1, format="csr", random_state=rng)

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # Add a layer
    adata.layers["log1p"] = X.copy()
    adata.layers["log1p"].data = np.log1p(adata.layers["log1p"].data)

    return adata
