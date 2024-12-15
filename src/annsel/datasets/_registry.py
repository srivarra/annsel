from enum import Enum

import pooch

ANNSEL_OS_CACHE = "annsel-cxg-datasets"

CZI_DATASETS_BASE_URL = "https://datasets.cellxgene.cziscience.com/"


class CXRDatasetsRegistry(str, Enum):
    """A collection of cell x gene datasets."""

    LEUKEMIC_BONE_MARROW_DONORS = "a48b7e8a-9db3-45f1-9729-c87738c0082f.h5ad"
    "15 leukemic bone marrow donors"

    FOCAL_CORTICAL_DYSPLASIA_DONORS = "29d50da1-430a-482c-ac9a-7077c727aa34.h5ad"
    "8 focal cortical dysplasia donors"

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
    base_url=CZI_DATASETS_BASE_URL,
    registry={
        CXRDatasetsRegistry.LEUKEMIC_BONE_MARROW_DONORS.fname: "e09ea1acda5c7ca39ae42da134ea8cf9e553d99931ee3d6372316ff0137dc93a",
        CXRDatasetsRegistry.FOCAL_CORTICAL_DYSPLASIA_DONORS.fname: "b22c63cb8665c77ad041f946c733ff9f7ff53008b4a0a59a2294a141bf31355f",
    },
)
