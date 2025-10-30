"""Tests for with_var() - gene-level cross-component operations."""

import anndata as ad
import narwhals as nw
import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def adata_gene_analysis():
    """Create AnnData for gene-level analysis."""
    # Cell types for representative samples
    obs = pd.DataFrame(
        {
            "cell_type": ["B_cell", "T_cell", "NK_cell", "Monocyte"],
            "sample": ["S1", "S1", "S2", "S2"],
        },
        index=["B_cell_1", "T_cell_1", "NK_cell_1", "Mono_1"],
    )

    # Gene metadata with biotypes
    var = pd.DataFrame(
        {
            "gene_type": ["protein_coding", "lncRNA", "protein_coding", "protein_coding", "miRNA"],
            "chromosome": ["chr1", "chr2", "chr1", "chrX", "chr3"],
            "highly_variable": [True, False, True, False, False],
        },
        index=["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
    )

    # Expression matrix (genes × cells perspective)
    X = np.array(
        [
            [10.0, 5.0, 2.0, 8.0],  # GENE1 - high in B and Mono
            [2.0, 3.0, 1.5, 2.5],  # GENE2 - lncRNA, low overall
            [15.0, 12.0, 18.0, 10.0],  # GENE3 - high variable
            [5.0, 5.0, 5.0, 5.0],  # GENE4 - constant expression
            [0.5, 0.8, 0.3, 0.6],  # GENE5 - miRNA, very low
        ]
    )

    return ad.AnnData(X=X.T, obs=obs, var=var)  # Note: X.T for correct orientation


def test_with_var_basic(adata_gene_analysis) -> None:
    """Test combining var with expression across cells."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_gene_analysis)

    # Combine gene metadata with expression in specific cells
    combined = nwdata.with_var(
        "gene_type",
        "highly_variable",
        obs_names=["B_cell_1", "T_cell_1"],
    )

    df = combined.to_native()

    # Verify structure
    assert "gene_type" in df.columns
    assert "highly_variable" in df.columns
    assert "B_cell_1" in df.columns
    assert "T_cell_1" in df.columns
    assert len(df) == 5  # 5 genes


def test_gene_type_expression_stats(adata_gene_analysis) -> None:
    """Test grouping genes by type and computing expression stats."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_gene_analysis)

    # Combine var with expression in all cells
    combined = nwdata.with_var(
        "gene_type",
        obs_names=["B_cell_1", "T_cell_1", "NK_cell_1", "Mono_1"],
    )

    # Group by gene type and compute stats
    stats = combined.group_by("gene_type").agg(
        [
            nw.col("B_cell_1").mean().alias("mean_in_B"),
            nw.col("T_cell_1").mean().alias("mean_in_T"),
            nw.col("NK_cell_1").mean().alias("mean_in_NK"),
            nw.col("Mono_1").mean().alias("mean_in_Mono"),
            nw.len().alias("n_genes"),
        ]
    )

    result = stats.to_native()

    # Verify grouping
    assert "lncRNA" in result["gene_type"].values
    assert "protein_coding" in result["gene_type"].values
    assert "miRNA" in result["gene_type"].values

    # Protein coding should have 3 genes
    pc_row = result[result["gene_type"] == "protein_coding"].iloc[0]
    assert pc_row["n_genes"] == 3


def test_chromosome_expression_analysis(adata_gene_analysis) -> None:
    """Test grouping genes by chromosome and computing expression."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_gene_analysis)

    # Combine var with expression
    combined = nwdata.with_var(
        "chromosome",
        obs_names=["B_cell_1", "T_cell_1"],
    )

    # Group by chromosome
    chr_stats = combined.group_by("chromosome").agg(
        [
            nw.col("B_cell_1").mean().alias("B_mean"),
            nw.col("T_cell_1").mean().alias("T_mean"),
            nw.len().alias("n_genes"),
        ]
    )

    result = chr_stats.to_native()

    # chr1 should have 2 genes (GENE1, GENE3)
    chr1 = result[result["chromosome"] == "chr1"].iloc[0]
    assert chr1["n_genes"] == 2


def test_highly_variable_expression_patterns(adata_gene_analysis) -> None:
    """Test comparing expression of highly variable vs other genes."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_gene_analysis)

    # Combine all var columns with all cell expression
    combined = nwdata.with_var("highly_variable")

    # Compare highly variable vs not
    stats = combined.group_by("highly_variable").agg(
        [
            nw.col("B_cell_1").mean().alias("B_mean"),
            nw.col("T_cell_1").mean().alias("T_mean"),
            nw.col("NK_cell_1").mean().alias("NK_mean"),
            nw.col("B_cell_1").std().alias("B_std"),
            nw.col("T_cell_1").std().alias("T_std"),
            nw.len().alias("n_genes"),
        ]
    )

    result = stats.to_native()

    # Should have 2 groups (True, False)
    assert len(result) == 2
    hv_row = result[result["highly_variable"] == True].iloc[0]  # noqa: E712
    assert hv_row["n_genes"] == 2  # GENE1 and GENE3


def test_with_var_all_columns(adata_gene_analysis) -> None:
    """Test with_var without specifying var_cols (includes all)."""
    from narwhals._utils import Version

    from annsel.narwhals_plugin import __narwhals_namespace__

    namespace = __narwhals_namespace__(Version.MAIN)
    nwdata = namespace.from_native(adata_gene_analysis)

    # Include all var columns
    combined = nwdata.with_var(obs_names=["B_cell_1", "T_cell_1"])

    df = combined.to_native()

    # Should have all var columns
    assert "gene_type" in df.columns
    assert "chromosome" in df.columns
    assert "highly_variable" in df.columns

    # Plus the selected cells
    assert "B_cell_1" in df.columns
    assert "T_cell_1" in df.columns


def test_with_var_integration_with_annsel(adata_gene_analysis) -> None:
    """Test with_var through annsel .an accessor."""
    # Use through annsel (unified API)
    combined = adata_gene_analysis.an.with_var(
        "gene_type",
        obs_names=["B_cell_1", "T_cell_1"],
    )

    # Group by gene type
    stats = combined.group_by("gene_type").agg(
        [
            nw.col("B_cell_1").mean().alias("B_mean"),
            nw.col("T_cell_1").mean().alias("T_mean"),
        ]
    )

    result = stats.to_native()
    assert len(result) == 3  # protein_coding, lncRNA, miRNA
