from functools import reduce
from operator import and_
from typing import Literal

import anndata as ad
import pandas as pd
from scipy import sparse


def _get_final_indices(
    obj_names: pd.Index,
    *idx: pd.Index,
) -> pd.Index:
    return obj_names.intersection(pd.Index(list(reduce(and_, map(set, *idx)))))


def _construct_adata_from_indices(
    adata: ad.AnnData,
    obs_idx: pd.Index,
    var_idx: pd.Index,
    sparse_method: Literal["csr", "csc"] | str | None = None,
) -> ad.AnnData:
    _adata = adata[obs_idx, var_idx]

    match sparse_method:
        case "csc":
            _X = sparse.csc_matrix(_adata.X)
        case "csr":
            _X = sparse.csr_matrix(_adata.X)
        case _:
            _X = _adata.X

    return _adata
