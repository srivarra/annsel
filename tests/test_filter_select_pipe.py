import anndata as ad
import anndata.tests.helpers as ath
import numpy as np
import pytest

import annsel as an


@pytest.mark.filter
class TestFilterAnnData:
    def test_filter_var_names(self, lbm_dataset: ad.AnnData) -> None:
        """Test filtering by variable names."""
        adata = lbm_dataset.an.filter(var_names=an.var_names.str.starts_with("ENSG0000018"))
        verify_adata = lbm_dataset[:, lbm_dataset.var_names.str.startswith("ENSG0000018")]

        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_var(self, lbm_dataset: ad.AnnData) -> None:
        """Test filtering by variable column."""
        adata = lbm_dataset.an.filter(var=an.col(["vst.mean"]) >= 1)
        verify_adata = lbm_dataset[:, lbm_dataset.var["vst.mean"] >= 1]

        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_obs(self, lbm_dataset: ad.AnnData) -> None:
        """Test filtering by observation column."""
        adata = lbm_dataset.an.filter(obs=an.col("Cell_label") == "Lymphomyeloid prog")
        verify_adata = lbm_dataset[lbm_dataset.obs["Cell_label"] == "Lymphomyeloid prog"]

        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_obs_names(self, lbm_dataset: ad.AnnData) -> None:
        """Test filtering by observation names."""
        adata = lbm_dataset.an.filter(obs_names=an.obs_names.str.ends_with("4"))
        verify_adata = lbm_dataset[lbm_dataset.obs_names.str.endswith("4")]

        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_X_sparse(self, lbm_dataset: ad.AnnData) -> None:
        """Test filtering X with sparse matrix (default case)."""
        adata = lbm_dataset.an.filter(x=an.col(["ENSG00000206560"]) >= 1)
        verify_adata = lbm_dataset[lbm_dataset[:, "ENSG00000206560"].X >= 1, :]
        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_X_dense(self, lbm_dataset: ad.AnnData) -> None:
        """Test filtering X with dense matrix."""
        # Convert to dense
        dense_adata = lbm_dataset.copy()
        dense_adata.X = dense_adata.X.toarray()

        adata = dense_adata.an.filter(x=an.col(["ENSG00000206560"]) >= 1)
        verify_adata = dense_adata[dense_adata[:, "ENSG00000206560"].X >= 1, :]
        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_X_layer(self, lbm_dataset: ad.AnnData) -> None:
        """Test filtering X from a specific layer."""
        # Add a test layer
        lbm_dataset.layers["test"] = np.log1p(lbm_dataset.X.toarray())

        adata = lbm_dataset.an.filter(x=an.col(["ENSG00000206560"]) >= 1, layer="test")
        verify_adata = lbm_dataset[lbm_dataset[:, "ENSG00000206560"].layers["test"] >= 1, :]
        ath.assert_adata_equal(adata, verify_adata)


@pytest.mark.pipe
class TestPipeAnnData:
    def test_pipe_simple(self, lbm_dataset: ad.AnnData) -> None:
        """Test basic pipe functionality."""

        def count_cells(adata):
            return adata.n_obs

        result = lbm_dataset.an.pipe(count_cells)
        assert result == lbm_dataset.n_obs

    def test_pipe_with_kwargs(self, lbm_dataset: ad.AnnData) -> None:
        """Test pipe with kwargs."""

        def subset_and_count(adata, cell_type):
            return adata[adata.obs["Cell_label"] == cell_type].n_obs

        result = lbm_dataset.an.pipe(subset_and_count, cell_type="Lymphomyeloid prog")
        verify_count = len(lbm_dataset[lbm_dataset.obs["Cell_label"] == "Lymphomyeloid prog"])
        assert result == verify_count

    def test_pipe_with_target(self, lbm_dataset: ad.AnnData) -> None:
        """Test pipe with explicit target argument."""

        def custom_filter(data_key, cell_type):
            return data_key[data_key.obs["Cell_label"] == cell_type]

        result = lbm_dataset.an.pipe((custom_filter, "data_key"), cell_type="Lymphomyeloid prog")

        verify_adata = lbm_dataset[lbm_dataset.obs["Cell_label"] == "Lymphomyeloid prog"]
        ath.assert_adata_equal(result, verify_adata)

    def test_pipe_chaining(self, lbm_dataset: ad.AnnData) -> None:
        """Test chaining multiple pipe operations."""

        def subset_by_cell_type(adata, cell_type):
            return adata[adata.obs["Cell_label"] == cell_type]

        def get_obs_names(adata):
            return adata.obs_names.tolist()

        result = lbm_dataset.an.pipe(subset_by_cell_type, cell_type="Lymphomyeloid prog").an.pipe(get_obs_names)

        verify_result = lbm_dataset[lbm_dataset.obs["Cell_label"] == "Lymphomyeloid prog"].obs_names.tolist()
        assert result == verify_result

    def test_pipe_raises_on_duplicate_target(self, lbm_dataset: ad.AnnData) -> None:
        """Test error handling for duplicate targets."""

        def dummy_func(adata_param):
            return adata_param

        with pytest.raises(ValueError, match="adata_param is both the pipe target and a keyword argument"):
            lbm_dataset.an.pipe((dummy_func, "adata_param"), adata_param=lbm_dataset)


@pytest.mark.select
class TestSelectAnnData:
    def test_select_with_single_predicate(self, lbm_dataset: ad.AnnData) -> None:
        """Test selecting observations with a single predicate."""
        adata = lbm_dataset.an.select(obs=an.col(["Cell_label", "is_primary_data"]))

        verify_adata = lbm_dataset.copy()
        verify_adata.obs = verify_adata.obs[["Cell_label", "is_primary_data"]]
        ath.assert_adata_equal(adata, verify_adata)

    def test_select_with_multiple_predicates(self, lbm_dataset: ad.AnnData) -> None:
        """Test selecting observations and variables with multiple predicates."""
        adata = lbm_dataset.an.select(obs=an.col(["Cell_label"]), var=an.col(["vst.mean", "feature_type"]))

        verify_adata = lbm_dataset.copy()
        verify_adata.obs = verify_adata.obs[["Cell_label"]]
        verify_adata.var = verify_adata.var[["vst.mean", "feature_type"]]
        ath.assert_adata_equal(adata, verify_adata, exact=False)

    def test_select_with_x_predicate(self, lbm_dataset: ad.AnnData) -> None:
        """Test selecting observations with an X predicate."""
        adata = lbm_dataset.an.select(x=an.col(["ENSG00000204472", "ENSG00000206560"]))

        verify_adata = lbm_dataset[:, ["ENSG00000204472", "ENSG00000206560"]]

        ath.assert_adata_equal(adata, verify_adata)
