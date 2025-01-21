from dataclasses import dataclass
from enum import Enum
from typing import Literal

import anndata as ad
import pooch
from upath import UPath

ANNSEL_OS_CACHE = "annsel-cxg-datasets"

CZI_DATASETS_BASE_URL = "https://datasets.cellxgene.cziscience.com/"


@dataclass
class CXRDataset:
    """A dataclass for a dataset from the CellXGene repository.

    Attributes
    ----------
    name
        The filename of the dataset.
    sha256
        The SHA256 checksum of the dataset.
    """

    name: str
    sha256: str


class CXRDatasetsRegistry(Enum):
    """A collection of cell x gene datasets."""

    LEUKEMIC_BONE_MARROW_DONORS = CXRDataset(
        name="a48b7e8a-9db3-45f1-9729-c87738c0082f.h5ad",
        sha256="e09ea1acda5c7ca39ae42da134ea8cf9e553d99931ee3d6372316ff0137dc93a",
    )
    """
    Leukemic Bone Marrow Donors dataset.
    """

    @property
    def fname(self) -> str:
        """Return the filename.

        Returns
        -------
        The record's filename.
        """
        return self.value.name

    @property
    def sha256(self) -> str:
        """Return the SHA256 checksum.

        Returns
        -------
        The record's SHA256 checksum.
        """
        return self.value.sha256

    def fetch(self) -> UPath:
        """Download and return the path to the dataset.

        Returns
        -------
        The absolute path to the dataset.
        """
        p = pooch.create(
            path=pooch.os_cache(ANNSEL_OS_CACHE),
            base_url=CZI_DATASETS_BASE_URL,
            registry={
                self.fname: self.sha256,
            },
        )
        return UPath(p.fetch(self.fname, progressbar=True))


def leukemic_bone_marrow_dataset(
    backed: bool | Literal["r", "r+"] | None = None,
) -> ad.AnnData:
    """Load a leukemic bone marrow cytometry dataset.

    This single cell cytometry dataset is from 15 leukemic bone marrow donors. It consists of
    31,586 cells and 458 genes.

    The saved file contains the following notable metrics (and more):
        - `obs`: clustering information, cell type associations, developmental stage,
        - `var`: vst.mean, vst.variance, feature_name, feature_type
        - `uns`: cell_type_ontology_term_id_colors
        - `obsm`: PCA, UMAP and tSNE embeddings

    Parameters
    ----------
    backed
        Whether to load the dataset in backed mode.

    Notes
    -----
    View the dataset at `CxG`_.

    References
    ----------
    :cite:p:`Triana2021`

    Download Size: 25.7 MB

    .. _CxG: https://cellxgene.cziscience.com/e/b3a5a10f-b1cb-4e8e-abce-bf345448625b.cxg/

    Returns
    -------
    The Leukemic Bone Marrow Donors dataset as an AnnData object.

    """
    _adata_path: UPath = CXRDatasetsRegistry.LEUKEMIC_BONE_MARROW_DONORS.fetch()

    adata = ad.io.read_h5ad(
        _adata_path,
        backed=backed,
    )
    adata.strings_to_categoricals()
    return adata


def marimo_dataset() -> ad.AnnData:
    """Generate a fake in memory dataset for Marimo, since it's not possible to load a dataset from a URL.

    Returns
    -------
    A minimal  dataset as an AnnData object.
    """
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(42)

    X = rng.random((100, 100)) * 10
    log1p_X = np.log1p(X)

    obs = pd.DataFrame(
        data={
            "cell_type": pd.Categorical(
                rng.choice(a=["A", "B", "C"], p=[0.3, 0.4, 0.3], size=100),
            ),
            "batch": pd.Categorical(
                rng.choice(a=["1", "2", "3"], p=[0.3, 0.4, 0.3], size=100),
            ),
            "development_stage": pd.Categorical(
                rng.choice(a=["embryonic", "adult"], p=[0.3, 0.7], size=100),
            ),
        },
        index=pd.Index([f"cell_{i}" for i in range(100)], name="index"),
    )

    var = pd.DataFrame(
        data={
            "feature_type": pd.Categorical(
                rng.choice(a=["protein", "rna"], p=[0.3, 0.7], size=100),
            ),
        },
        index=pd.Index([f"gene_{i}" for i in range(100)], name="index"),
    )

    adata = ad.AnnData(X=X, obs=obs, var=var, layers={"log1p": log1p_X})
    return adata

    # adata = ad.AnnData(X=X, obs=obs, var=var)
    # return adata
