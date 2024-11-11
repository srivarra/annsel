import anndata as ad
import anndata.tests.helpers as ath
import pytest
from scipy import sparse

import annsel as an


class TestFilterAnnData:
    def test_filter_var_names(self, lbm_dataset: ad.AnnData):
        adata = lbm_dataset.an.filter(an.var_names().str.starts_with("ENSG0000018"))
        verify_adata = lbm_dataset[:, lbm_dataset.var_names.str.startswith("ENSG0000018")]

        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_var(self, lbm_dataset: ad.AnnData):
        adata = lbm_dataset.an.filter(an.var_col(["vst.mean"]) >= 1)
        verify_adata = lbm_dataset[:, lbm_dataset.var["vst.mean"] >= 1]

        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_obs(self, lbm_dataset: ad.AnnData):
        adata = lbm_dataset.an.filter(an.obs_col("Cell_label") == "Lymphomyeloid prog")
        verify_adata = lbm_dataset[lbm_dataset.obs["Cell_label"] == "Lymphomyeloid prog"]

        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_obs_names(self, lbm_dataset: ad.AnnData):
        adata = lbm_dataset.an.filter(an.obs_names().str.ends_with("4"))
        verify_adata = lbm_dataset[lbm_dataset.obs_names.str.endswith("4")]

        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_X_sparse(self, lbm_dataset: ad.AnnData):
        """Test filtering X with sparse matrix (default case)."""
        adata = lbm_dataset.an.filter(an.x(["ENSG00000206560"]) >= 1)
        verify_adata = lbm_dataset[lbm_dataset[:, "ENSG00000206560"].X >= 1, :]
        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_X_dense(self, lbm_dataset: ad.AnnData):
        """Test filtering X with dense matrix."""
        # Convert to dense
        dense_adata = lbm_dataset.copy()
        dense_adata.X = dense_adata.X.toarray()

        adata = dense_adata.an.filter(an.x(["ENSG00000206560"]) >= 1)
        verify_adata = dense_adata[dense_adata[:, "ENSG00000206560"].X >= 1, :]
        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_X_force_csr(self, lbm_dataset: ad.AnnData):
        """Test filtering X with forced CSR format."""
        dense_adata = lbm_dataset.copy()
        dense_adata.X = dense_adata.X.toarray()

        adata = lbm_dataset.an.filter(an.x(["ENSG00000206560"]) >= 1, keep_sparse=True, sparse_method="csr")
        verify_adata = dense_adata[lbm_dataset[:, "ENSG00000206560"].X >= 1, :]
        ath.assert_adata_equal(adata, verify_adata)
        assert sparse.isspmatrix_csr(adata.X)

    def test_filter_X_force_csc(self, lbm_dataset: ad.AnnData):
        """Test filtering X with forced CSC format."""
        # Convert to dense first to test conversion
        dense_adata = lbm_dataset.copy()
        dense_adata.X = dense_adata.X.toarray()

        adata = dense_adata.an.filter(an.x(["ENSG00000206560"]) >= 1, keep_sparse=True, sparse_method="csc")
        verify_adata = dense_adata[dense_adata[:, "ENSG00000206560"].X >= 1, :]
        ath.assert_adata_equal(adata, verify_adata)
        assert sparse.isspmatrix_csc(adata.X)

    def test_filter_X_no_sparsity(self, lbm_dataset: ad.AnnData):
        """Test filtering X with sparsity disabled."""
        adata = lbm_dataset.an.filter(an.x(["ENSG00000206560"]) >= 1, keep_sparse=False)
        verify_adata = lbm_dataset[lbm_dataset[:, "ENSG00000206560"].X >= 1, :]
        ath.assert_adata_equal(adata, verify_adata)
        assert not sparse.issparse(adata.X)

    def test_filter_X_layer(self, lbm_dataset: ad.AnnData):
        """Test filtering X from a specific layer."""
        # Add a test layer
        lbm_dataset.layers["test"] = lbm_dataset.X.copy()

        adata = lbm_dataset.an.filter(an.x(["ENSG00000206560"]) >= 1, layer="test")
        verify_adata = lbm_dataset[lbm_dataset[:, "ENSG00000206560"].layers["test"] >= 1, :]
        ath.assert_adata_equal(adata, verify_adata)


class TestPipeAnnData:
    def test_pipe_simple(self, lbm_dataset: ad.AnnData):
        """Test basic pipe functionality."""

        def count_cells(adata):
            return adata.n_obs

        result = lbm_dataset.an.pipe(count_cells)
        assert result == lbm_dataset.n_obs

    def test_pipe_with_kwargs(self, lbm_dataset: ad.AnnData):
        """Test pipe with kwargs."""

        def subset_and_count(adata, cell_type):
            return adata[adata.obs["Cell_label"] == cell_type].n_obs

        result = lbm_dataset.an.pipe(subset_and_count, cell_type="Lymphomyeloid prog")
        verify_count = len(lbm_dataset[lbm_dataset.obs["Cell_label"] == "Lymphomyeloid prog"])
        assert result == verify_count

    def test_pipe_with_target(self, lbm_dataset: ad.AnnData):
        """Test pipe with explicit target argument."""

        def custom_filter(data_key, cell_type):
            return data_key[data_key.obs["Cell_label"] == cell_type]

        result = lbm_dataset.an.pipe((custom_filter, "data_key"), cell_type="Lymphomyeloid prog")

        verify_adata = lbm_dataset[lbm_dataset.obs["Cell_label"] == "Lymphomyeloid prog"]
        ath.assert_adata_equal(result, verify_adata)

    def test_pipe_chaining(self, lbm_dataset: ad.AnnData):
        """Test chaining multiple pipe operations."""

        def subset_by_cell_type(adata, cell_type):
            return adata[adata.obs["Cell_label"] == cell_type]

        def get_obs_names(adata):
            return adata.obs_names.tolist()

        result = lbm_dataset.an.pipe(subset_by_cell_type, cell_type="Lymphomyeloid prog").an.pipe(get_obs_names)

        verify_result = lbm_dataset[lbm_dataset.obs["Cell_label"] == "Lymphomyeloid prog"].obs_names.tolist()
        assert result == verify_result

    def test_pipe_raises_on_duplicate_target(self, lbm_dataset: ad.AnnData):
        """Test error handling for duplicate targets."""

        def dummy_func(adata_param):
            return adata_param

        with pytest.raises(ValueError, match="adata_param is both the pipe target and a keyword argument"):
            lbm_dataset.an.pipe((dummy_func, "adata_param"), adata_param=lbm_dataset)
