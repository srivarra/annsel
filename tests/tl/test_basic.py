import anndata as ad
import anndata.tests.helpers as ath
import pytest

import annsel as an


class TestFilter:
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

    @pytest.mark.xfail(reason="Issues with spare X")
    def test_filter_X(self, lbm_dataset: ad.AnnData):
        adata = lbm_dataset.an.filter(an.x(["HES4"]) >= 1)
        verify_adata = lbm_dataset[lbm_dataset[:, "HES4"].X >= 1, :]

        ath.assert_adata_equal(adata, verify_adata)
