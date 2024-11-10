import anndata as ad
import pytest

from annsel.datasets import leukemic_bone_marrow_dataset


@pytest.fixture
def lbm_dataset() -> ad.AnnData:
    return leukemic_bone_marrow_dataset()
