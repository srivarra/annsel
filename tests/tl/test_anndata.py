import anndata as ad
import anndata.tests.helpers as ath
import numpy as np
import pytest

import annsel as an


@pytest.mark.filter
class TestFilterAnnData:
    def test_filter_var_names(self, lbm_dataset: ad.AnnData):
        """Test filtering by variable names."""
        adata = lbm_dataset.an.filter(var_names=an.var_names.str.starts_with("ENSG0000018"))
        verify_adata = lbm_dataset[:, lbm_dataset.var_names.str.startswith("ENSG0000018")]

        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_var(self, lbm_dataset: ad.AnnData):
        """Test filtering by variable column."""
        adata = lbm_dataset.an.filter(var=an.col(["vst.mean"]) >= 1)
        verify_adata = lbm_dataset[:, lbm_dataset.var["vst.mean"] >= 1]

        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_obs(self, lbm_dataset: ad.AnnData):
        """Test filtering by observation column."""
        adata = lbm_dataset.an.filter(obs=an.col("Cell_label") == "Lymphomyeloid prog")
        verify_adata = lbm_dataset[lbm_dataset.obs["Cell_label"] == "Lymphomyeloid prog"]

        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_obs_names(self, lbm_dataset: ad.AnnData):
        """Test filtering by observation names."""
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

    def test_filter_X_layer(self, lbm_dataset: ad.AnnData):
        """Test filtering X from a specific layer."""
        # Add a test layer
        lbm_dataset.layers["test"] = np.log1p(lbm_dataset.X.toarray())

        adata = lbm_dataset.an.filter(x=an.col(["ENSG00000206560"]) >= 1, layer="test")
        verify_adata = lbm_dataset[lbm_dataset[:, "ENSG00000206560"].layers["test"] >= 1, :]
        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_obsm(self, lbm_dataset: ad.AnnData, rng: np.random.Generator):
        """Test filtering based on obsm matrices."""
        # Create a test obsm key if it doesn't exist
        if "X_test" not in lbm_dataset.obsm:
            lbm_dataset.obsm["X_test"] = rng.random((lbm_dataset.n_obs, 2))
            # Make some values specifically >0.8 for testing
            lbm_dataset.obsm["X_test"][0:5, 0] = 0.9

        # Filter cells where first coordinate in X_test > 0.8
        adata = lbm_dataset.an.filter(obsm={"X_test": an.col([0]) > 0.8})

        # Verify the filter worked correctly
        cells_with_high_values = lbm_dataset.obs_names[lbm_dataset.obsm["X_test"][:, 0] > 0.8]
        verify_adata = lbm_dataset[cells_with_high_values]

        ath.assert_adata_equal(adata, verify_adata)

    def test_filter_varm(self, lbm_dataset: ad.AnnData, rng: np.random.Generator):
        """Test filtering based on varm matrices."""
        # Create a test varm key if it doesn't exist
        if "PCs" not in lbm_dataset.varm:
            lbm_dataset.varm["PCs"] = rng.random((lbm_dataset.n_vars, 3))
            # Make some values specifically <0.2 for testing
            lbm_dataset.varm["PCs"][0:10, 0] = 0.1

        # Filter genes where first PC coordinate in PCs < 0.2
        adata = lbm_dataset.an.filter(varm={"PCs": an.col([0]) < 0.2})

        # Verify the filter worked correctly
        genes_with_low_values = lbm_dataset.var_names[lbm_dataset.varm["PCs"][:, 0] < 0.2]
        verify_adata = lbm_dataset[:, genes_with_low_values]

        ath.assert_adata_equal(adata, verify_adata)


@pytest.mark.pipe
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


@pytest.mark.select
class TestSelectAnnData:
    def test_select_with_single_predicate(self, lbm_dataset: ad.AnnData):
        """Test selecting observations with a single predicate."""
        adata = lbm_dataset.an.select(obs=an.col(["Cell_label", "is_primary_data"]))

        verify_adata = lbm_dataset.copy()
        verify_adata.obs = verify_adata.obs[["Cell_label", "is_primary_data"]]
        ath.assert_adata_equal(adata, verify_adata)

    def test_select_with_multiple_predicates(self, lbm_dataset: ad.AnnData):
        """Test selecting observations and variables with multiple predicates."""
        adata = lbm_dataset.an.select(obs=an.col(["Cell_label"]), var=an.col(["vst.mean", "feature_type"]))

        verify_adata = lbm_dataset.copy()
        verify_adata.obs = verify_adata.obs[["Cell_label"]]
        verify_adata.var = verify_adata.var[["vst.mean", "feature_type"]]
        ath.assert_adata_equal(adata, verify_adata, exact=False)

    def test_select_with_x_predicate(self, lbm_dataset: ad.AnnData):
        """Test selecting observations with an X predicate."""
        adata = lbm_dataset.an.select(x=an.col(["ENSG00000204472", "ENSG00000206560"]))

        verify_adata = lbm_dataset[:, ["ENSG00000204472", "ENSG00000206560"]]

        ath.assert_adata_equal(adata, verify_adata)


@pytest.mark.groupby
class TestGroupByAnnData:
    def test_group_by_obs(self, lbm_dataset: ad.AnnData):
        """Test grouping by observation columns."""
        groups = list(lbm_dataset.an.group_by(obs=an.col("Cell_label")))

        # Verify number of groups matches unique Cell_label values
        assert len(groups) == lbm_dataset.obs["Cell_label"].nunique()

        # Verify each group is a GroupByAnndata instance
        for group in groups:
            # Verify group contains correct data
            cell_type = group.adata.obs["Cell_label"].iloc[0]
            verify_adata = lbm_dataset[lbm_dataset.obs["Cell_label"] == cell_type]
            ath.assert_adata_equal(group.adata, verify_adata)

            # Verify obs_dict contains the expected values
            assert group.obs_dict == {"Cell_label": cell_type}

    def test_group_by_var(self, lbm_dataset: ad.AnnData):
        """Test grouping by variable columns."""
        groups = list(lbm_dataset.an.group_by(var=an.col("feature_type")))

        # Verify number of groups matches unique feature_type values
        assert len(groups) == lbm_dataset.var["feature_type"].nunique()

        # Verify each group is a GroupByAnndata instance
        for group in groups:
            # Verify group contains correct data
            feature_type = group.adata.var["feature_type"].iloc[0]
            verify_adata = lbm_dataset[:, lbm_dataset.var["feature_type"] == feature_type]
            ath.assert_adata_equal(group.adata, verify_adata)

            # Verify var_dict contains the expected values
            assert group.var_dict == {"feature_type": feature_type}

    def test_group_by_multiple_columns(self, lbm_dataset: ad.AnnData):
        """Test grouping by multiple columns simultaneously."""
        groups = list(lbm_dataset.an.group_by(obs=an.col(["Cell_label", "sex"]), copy=True))

        # Verify number of groups matches unique combinations
        expected_groups = lbm_dataset.obs.groupby(["Cell_label", "sex"], observed=True).ngroups
        assert len(groups) == expected_groups

        # Verify each group has the correct columns in obs_dict
        for group in groups:
            assert set(group.obs_dict.keys()) == {"Cell_label", "sex"}

    def test_group_by_no_kwargs(self, lbm_dataset: ad.AnnData):
        """Test grouping with no kwargs returns the original AnnData."""
        # When no grouping arguments are given, it should return the original adata
        # (or a copy if copy=True)
        grouped_obj = lbm_dataset.an.group_by()
        assert isinstance(grouped_obj, ad.AnnData)
        ath.assert_adata_equal(grouped_obj, lbm_dataset)

        # With copy=True
        copied_grouped_obj = lbm_dataset.an.group_by(copy=True)
        assert isinstance(copied_grouped_obj, ad.AnnData)
        ath.assert_adata_equal(copied_grouped_obj, lbm_dataset)
        # Ensure it's a copy, not the same object
        assert copied_grouped_obj is not lbm_dataset

    def test_group_by_obs_and_var(self, lbm_dataset: ad.AnnData):
        """Test grouping by both observations and variables."""
        groups = list(lbm_dataset.an.group_by(obs=an.col(["Cell_label"]), var=an.col(["feature_type"])))

        # Verify structure with both obs and var grouping
        for group in groups:
            # Check that obs and var dictionaries contain the expected keys
            assert "Cell_label" in group.obs_dict
            assert "feature_type" in group.var_dict

            # Verify the subset contains the expected data
            cell_type = group.obs_dict["Cell_label"]
            feature_type = group.var_dict["feature_type"]

            verify_adata = lbm_dataset[
                lbm_dataset.obs["Cell_label"] == cell_type, lbm_dataset.var["feature_type"] == feature_type
            ]

            ath.assert_adata_equal(group.adata, verify_adata)

    def test_repr_format(self, lbm_dataset: ad.AnnData):
        """Test the tree-like representation of GroupByAnndata."""
        # Test with both obs and var specified
        group = next(lbm_dataset.an.group_by(obs=an.col(["Cell_label"]), var=an.col(["feature_type"])))

        # Check that the representation is formatted correctly
        repr_str = repr(group)

        # Check structure elements
        assert repr_str.startswith("GroupByAnnData:")
        assert "├── Observations:" in repr_str
        assert "├── Variables:" in repr_str
        assert "└── AnnData:" in repr_str

        # Check if content is included
        cell_type = group.obs_dict["Cell_label"]
        feature_type = group.var_dict["feature_type"]

        assert f"Cell_label: {cell_type}" in repr_str
        assert f"feature_type: {feature_type}" in repr_str

    def test_repr_format_no_obs(self, lbm_dataset: ad.AnnData):
        """Test representation when obs is None (all observations)."""
        # Test with only var specified (obs=None)
        group = next(lbm_dataset.an.group_by(var=an.col(["feature_type"])))

        repr_str = repr(group)

        # Check structure elements
        assert repr_str.startswith("GroupByAnnData:")
        assert "├── Observations:" in repr_str
        assert "(all observations)" in repr_str
        assert "├── Variables:" in repr_str

        # Check that var content is included
        feature_type = group.var_dict["feature_type"]
        assert f"feature_type: {feature_type}" in repr_str

    def test_repr_format_no_var(self, lbm_dataset: ad.AnnData):
        """Test representation when var is None (all variables)."""
        # Test with only obs specified (var=None)
        group = next(lbm_dataset.an.group_by(obs=an.col(["Cell_label"])))

        repr_str = repr(group)

        # Check structure elements
        assert repr_str.startswith("GroupByAnnData:")
        assert "├── Observations:" in repr_str
        assert "├── Variables:" in repr_str
        assert "(all variables)" in repr_str

        # Check that obs content is included
        cell_type = group.obs_dict["Cell_label"]
        assert f"Cell_label: {cell_type}" in repr_str

    def test_repr_format_multiple_groups(self, lbm_dataset: ad.AnnData):
        """Test representation when multiple groups are present."""
        # Create a new AnnData with multiple groups
        group = next(
            lbm_dataset.an.group_by(
                obs=an.col(["Cell_label", "disease"]), var=an.col(["feature_type", "feature_is_filtered"])
            )
        )

        repr_str = repr(group)

        assert repr_str.startswith("GroupByAnnData:")
        assert "├── Observations:" in repr_str
        assert "├── Variables:" in repr_str
        assert "└── AnnData:" in repr_str

        cell_type = group.obs_dict["Cell_label"]
        disease = group.obs_dict["disease"]
        feature_type = group.var_dict["feature_type"]
        feature_is_filtered = group.var_dict["feature_is_filtered"]

        assert f"Cell_label: {cell_type}" in repr_str
        assert f"disease: {disease}" in repr_str
        assert f"feature_type: {feature_type}" in repr_str
        assert f"feature_is_filtered: {feature_is_filtered}" in repr_str
