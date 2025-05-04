import anndata as ad
import narwhals as nw
import pandas as pd
import pytest

from annsel.core.utils import _construct_adata_from_indices, _extract_names_from_expr, _get_final_indices, second


def test_extract_names_from_expr():
    expr = nw.col(["a", "b"])
    assert _extract_names_from_expr(expr) == ("a", "b")

    assert _extract_names_from_expr("a") == ("a",)

    with pytest.raises(ValueError, match="Unknown predicate type"):
        _extract_names_from_expr(1)  # type: ignore


def test_get_final_indices():
    indices = [pd.Index([1, 2, 3]), pd.Index([1, 2, 6])]
    assert _get_final_indices(pd.Index([1, 2, 3, 4, 5, 6]), *indices).equals(pd.Index([1, 2]))

    assert _get_final_indices(pd.Index([1]), []).equals(pd.Index([]))


def test_second():
    assert second([1, 2, 3]) == 2
    assert second((1, 2, 3)) == 2
    assert second(pd.Index([1, 2, 3])) == 2
    assert second(pd.Index([1, 2, 3])) == 2

    with pytest.raises(ValueError, match="No second element found"):
        second([1])


def test_construct_adata_from_indices(rng):
    adata = ad.AnnData(X=rng.normal(size=(10, 10)))
    obs_idx = pd.Index([1, 2, 3])
    var_idx = pd.Index([4, 5, 6])
    adata_subset = _construct_adata_from_indices(adata, obs_idx, var_idx)
    assert adata_subset.n_obs == 3
    assert adata_subset.n_vars == 3
