"""Select operations using narwhals-anndata plugin backend."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import anndata as ad

    from annsel.core.typing import Predicates


def _select_plugin(
    adata: ad.AnnData,
    obs: Predicates | None = None,
    var: Predicates | None = None,
    x: Predicates | None = None,
    layer: str | None = None,
) -> ad.AnnData:
    """Select columns using the narwhals-anndata plugin.

    Parameters
    ----------
    adata
        The AnnData object.
    obs
        Column expressions to select from obs.
    var
        Column expressions to select from var.
    x
        Column expressions to select from X (by gene names).

    Returns
    -------
    ad.AnnData
        AnnData with selected columns.
    """
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    # Get plugin namespace with current narwhals version
    namespace = __narwhals_namespace__(Version.MAIN)
    result_adata = adata

    # Select obs columns
    if obs is not None:
        nwdata = namespace.from_native(result_adata)
        result_adata = nwdata.obs.select(obs)

    # Select var columns
    if var is not None:
        nwdata = namespace.from_native(result_adata)
        result_adata = nwdata.var.select(var)

    # Select columns from X (variables by var_names/index)
    if x is not None:
        # X selection means selecting specific genes by their var_names
        # Extract gene names and slice AnnData
        from annsel.core.utils import _extract_names_from_expr

        gene_names = _extract_names_from_expr(x)
        result_adata = result_adata[:, list(gene_names)]

    return result_adata
