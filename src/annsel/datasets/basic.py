from enum import Enum

import anndata as ad
import pooch
from upath import UPath

ANNSEL_OS_CACHE = "annsel-cxg-datasets"


class CXRDatasetsRegistry(str, Enum):
    """A collection of cell x gene datasets."""

    LEUKEMIC_BONE_MARROW_DONORS = "a48b7e8a-9db3-45f1-9729-c87738c0082f.h5ad"
    "15 leukemic bone marrow donors"

    @property
    def fname(self) -> str:
        """Return the filename.

        Returns
        -------
        The record's filename.
        """
        return self.value


cxg_pooch = pooch.create(
    path=pooch.os_cache(ANNSEL_OS_CACHE),
    base_url=UPath("https://datasets.cellxgene.cziscience.com/").as_uri(),
    registry={
        CXRDatasetsRegistry.LEUKEMIC_BONE_MARROW_DONORS.fname: "e09ea1acda5c7ca39ae42da134ea8cf9e553d99931ee3d6372316ff0137dc93a"
    },
)


def leukemic_bone_marrow_dataset() -> ad.AnnData:
    """Load a leukemic bone marrow cytometry dataset.

    This single cell cytometry dataset is from 15 leukemic bone marrow donors. It consists of
    31,586 cells and 458 genes.

    The saved file contains the following notable metrics (and more):
        - `obs`: clustering information, cell type associations, developmental stage,
        - `var`: vst.mean, vst.variance, feature_name, feature_type
        - `uns`: cell_type_ontology_term_id_colors
        - `obsm`: PCA, UMAP and tSNE embeddings

    Notes
    -----
    View the dataset at `CxG`_.

    Citation: :cite:p:`triana_single-cell_2021`

    Download Size: 25.7 MB

    .. _CxG: https://cellxgene.cziscience.com/e/b3a5a10f-b1cb-4e8e-abce-bf345448625b.cxg/

    Returns
    -------
    The Leukemic Bone Marrow Donors dataset as an AnnData object.

    """
    _adata_path = UPath(cxg_pooch.fetch(CXRDatasetsRegistry.LEUKEMIC_BONE_MARROW_DONORS.fname, progressbar=True))

    adata = ad.io.read_h5ad(
        _adata_path,
    )
    adata.strings_to_categoricals()
    return adata
