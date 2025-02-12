from functools import reduce
from operator import and_

import anndata as ad
import pandas as pd


def _get_final_indices(
    obj_names: pd.Index,
    *idx: pd.Index,
) -> pd.Index:
    return obj_names.intersection(pd.Index(list(reduce(and_, map(set, *idx)))))


def _construct_adata_from_indices(
    adata: ad.AnnData,
    obs_idx: pd.Index,
    var_idx: pd.Index,
) -> ad.AnnData:
    _adata = adata[obs_idx, var_idx]

    return _adata
