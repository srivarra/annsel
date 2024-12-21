import anndata as ad
import anndata.tests.helpers as ath
import numpy as np
import pytest
from scipy import sparse

import annsel as an


class TestFilterAnnData:
    def test_filter_var_names(self, lbm_dataset: ad.AnnData):
        adata = lbm_dataset.an.filter(var_names=an.var_names.str.starts_with("ENSG0000018"))
        verify_adata = lbm_dataset[:, lbm_dataset.var_names.str.startswith("ENSG0000018")]

        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_var(self, lbm_dataset: ad.AnnData):
        adata = lbm_dataset.an.filter(var=an.col(["vst.mean"]) >= 1)
        verify_adata = lbm_dataset[:, lbm_dataset.var["vst.mean"] >= 1]

        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_obs(self, lbm_dataset: ad.AnnData):
        adata = lbm_dataset.an.filter(obs=an.col("Cell_label") == "Lymphomyeloid prog")
        verify_adata = lbm_dataset[lbm_dataset.obs["Cell_label"] == "Lymphomyeloid prog"]

        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_obs_names(self, lbm_dataset: ad.AnnData):
        adata = lbm_dataset.an.filter(obs_names=an.obs_names.str.ends_with("4"))
        verify_adata = lbm_dataset[lbm_dataset.obs_names.str.endswith("4")]

        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_X_sparse(self, lbm_dataset: ad.AnnData):
        """Test filtering X with sparse matrix (default case)."""
        adata = lbm_dataset.an.filter(x=an.col(["ENSG00000206560"]) >= 1)
        verify_adata = lbm_dataset[lbm_dataset[:, "ENSG00000206560"].X >= 1, :]
        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_X_dense(self, lbm_dataset: ad.AnnData):
        """Test filtering X with dense matrix."""
        # Convert to dense
        dense_adata = lbm_dataset.copy()
        dense_adata.X = dense_adata.X.toarray()

        adata = dense_adata.an.filter(x=an.col(["ENSG00000206560"]) >= 1)
        verify_adata = dense_adata[dense_adata[:, "ENSG00000206560"].X >= 1, :]
        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_X_force_csr(self, lbm_dataset: ad.AnnData):
        """Test filtering X with forced CSR format."""
        dense_adata = lbm_dataset.copy()
        dense_adata.X = dense_adata.X.toarray()

        adata = lbm_dataset.an.filter(x=an.col(["ENSG00000206560"]) >= 1, sparse="csr")
        verify_adata = dense_adata[lbm_dataset[:, "ENSG00000206560"].X >= 1, :]
        ath.assert_adata_equal(adata, verify_adata)
        assert sparse.isspmatrix_csr(adata.X)

    def test_filter_X_layer(self, lbm_dataset: ad.AnnData):
        """Test filtering X from a specific layer."""
        # Add a test layer
        lbm_dataset.layers["test"] = np.log(lbm_dataset.X.toarray())

        adata = lbm_dataset.an.filter(x=an.col(["ENSG00000206560"]) >= 1, layer="test")
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


class TestSelectAnnData:
    def test_select_with_single_predicate(self, lbm_dataset: ad.AnnData):
        adata = lbm_dataset.an.select(obs=an.col(["Cell_label", "is_primary_data"]))

        verify_adata = lbm_dataset.copy()
        verify_adata.obs = verify_adata.obs[["Cell_label", "is_primary_data"]]
        ath.assert_adata_equal(adata, verify_adata)

    def test_select_with_multiple_predicates(self, lbm_dataset: ad.AnnData):
        adata = lbm_dataset.an.select(obs=an.col(["Cell_label"]), var=an.col(["vst.mean", "feature_type"]))

        verify_adata = lbm_dataset.copy()
        verify_adata.obs = verify_adata.obs[["Cell_label"]]
        verify_adata.var = verify_adata.var[["vst.mean", "feature_type"]]
        ath.assert_adata_equal(adata, verify_adata, exact=False)

    def test_select_with_x_predicate(self, lbm_dataset: ad.AnnData):
        adata = lbm_dataset.an.select(x=an.col(["ENSG00000204472", "ENSG00000206560"]))

        verify_adata = lbm_dataset[:, ["ENSG00000204472", "ENSG00000206560"]]

        ath.assert_adata_equal(adata, verify_adata)


class TestGroupByAnnData:
    def test_group_by_obs(self, lbm_dataset: ad.AnnData):
        """Test grouping by observation columns."""
        groups = list(lbm_dataset.an.group_by(obs=an.col(["Cell_label"])))

        # Verify number of groups matches unique Cell_label values
        assert len(groups) == lbm_dataset.obs["Cell_label"].nunique()

        # Verify each group contains correct cells
        for group in groups:
            cell_type = group.obs["Cell_label"].iloc[0]
            verify_adata = lbm_dataset[lbm_dataset.obs["Cell_label"] == cell_type]
            ath.assert_adata_equal(group, verify_adata)

    def test_group_by_var(self, lbm_dataset: ad.AnnData):
        """Test grouping by variable columns."""
        groups = list(lbm_dataset.an.group_by(var=an.col(["feature_type"])))

        # Verify number of groups matches unique feature_type values
        assert len(groups) == lbm_dataset.var["feature_type"].nunique()

        # Verify each group contains correct features
        for group in groups:
            feature_type = group.var["feature_type"].iloc[0]
            verify_adata = lbm_dataset[:, lbm_dataset.var["feature_type"] == feature_type]
            ath.assert_adata_equal(group, verify_adata)

    def test_group_by_multiple_columns(self, lbm_dataset: ad.AnnData):
        """Test grouping by multiple columns simultaneously."""
        groups = list(lbm_dataset.an.group_by(obs=an.col(["Cell_label", "sex"])))

        # Verify number of groups matches unique combinations
        expected_groups = lbm_dataset.obs.groupby(["Cell_label", "sex"]).ngroups
        assert len(groups) == expected_groups

    def test_group_by_with_names(self, lbm_dataset: ad.AnnData):
        """Test grouping with return_group_names=True."""
        groups = list(lbm_dataset.an.group_by(obs=an.col(["Cell_label"]), return_group_names=True))

        # Verify structure of returned tuples
        for group_name, adata in groups:
            assert isinstance(group_name, tuple)
            assert isinstance(adata, ad.AnnData)
            # Verify group name matches the data
            assert adata.obs["Cell_label"].iloc[0] == group_name[0]

    def test_group_by_obs_and_var(self, lbm_dataset: ad.AnnData):
        """Test grouping by both observations and variables."""
        groups = list(
            lbm_dataset.an.group_by(obs=an.col(["Cell_label"]), var=an.col(["feature_type"]), return_group_names=True)
        )

        # Verify structure with both obs and var grouping
        for obs_name, var_name, adata in groups:
            assert isinstance(obs_name, tuple)
            assert isinstance(var_name, tuple)
            assert isinstance(adata, ad.AnnData)
            # Verify group names match the data
            assert adata.obs["Cell_label"].iloc[0] == obs_name[0]
            assert adata.var["feature_type"].iloc[0] == var_name[0]

    def test_group_by_empty_result(self, lbm_dataset: ad.AnnData):
        """Test grouping that results in no matches."""
        with pytest.raises(ValueError, match="No group keys passed!"):
            list(lbm_dataset.an.group_by(obs=an.col(["Cell_label"]) == "NonexistentType"))
