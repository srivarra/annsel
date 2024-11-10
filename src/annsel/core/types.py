from dataclasses import dataclass

import pandas as pd


@dataclass
class XIndicies:
    """A dataclass representing the indices of the X matrix of an AnnData object."""

    obs: pd.Index
    var: pd.Index
