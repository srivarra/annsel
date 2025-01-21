from annsel.datasets import marimo_dataset


def test_marimo_dataset():
    adata = marimo_dataset()
    assert adata.n_obs == 100
    assert adata.n_vars == 100
