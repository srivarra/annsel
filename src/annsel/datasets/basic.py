from typing import Literal

import anndata as ad
from upath import UPath

from ._registry import CXRDatasetsRegistry, cxg_pooch


def leukemic_bone_marrow_dataset(backed: bool | Literal["r", "r+"] | None = None) -> ad.AnnData:
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

    Citation: :cite:p:`Triana2021`

    Download Size: 25.7 MB

    .. _CxG: https://cellxgene.cziscience.com/e/b3a5a10f-b1cb-4e8e-abce-bf345448625b.cxg/

    Returns
    -------
    The Leukemic Bone Marrow Donors dataset as an AnnData object.

    """
    _adata_path: UPath = UPath(cxg_pooch.fetch(CXRDatasetsRegistry.LEUKEMIC_BONE_MARROW_DONORS.fname, progressbar=True))

    adata = ad.io.read_h5ad(
        _adata_path,
        backed=backed,
    )
    adata.strings_to_categoricals()
    return adata


def focal_cortical_dysplasia_dataset(backed: bool | Literal["r", "r+"] | None = None) -> ad.AnnData:
    """Load a dataset of focal cortical dysplasia donors.

    This snRNA-seq dataset is from 8 focal cortical dysplasia donors. It consists of
    61,525 cells and 36,406 genes.

    The saved file contains the following notable metrics (and more):
        - `obs`: donor_id, disease, developmental stage,
        - `var`: vst.mean, vst.variance, feature_name, feature_type
        - `uns`: cell_type_ontology_term_id_colors
        - `obsm`: PCA and UMAP embeddings

    Parameters
    ----------
    backed
        Whether to load the dataset in backed mode.

    Notes
    -----
    View the dataset at `CxG`_.

    Citation: :cite:p:`Galvo2024`

    Download Size: 681 MB

    .. _CxG: https://cellxgene.cziscience.com/e/01b84709-139c-4485-98a9-6e14c58fbbf6.cxg/

    Returns
    -------
    The Focal Cortical Dysplasia Donors dataset as an AnnData object.
    """
    _adata_path: UPath = UPath(
        cxg_pooch.fetch(CXRDatasetsRegistry.FOCAL_CORTICAL_DYSPLASIA_DONORS.fname, progressbar=True)
    )

    adata = ad.io.read_h5ad(
        _adata_path,
        backed=backed,
    )
    adata.strings_to_categoricals()
    return adata
