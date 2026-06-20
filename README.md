# annsel

<div align="center">

|               |                                                                                                                                                                                                              |
| :-----------: | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
|  **Status**   | [![Build][badge-build]][link-build] [![Tests][badge-test]][link-test] [![Documentation][badge-docs]][link-docs] [![codecov][badge-codecov]][link-codecov] [![prek][badge-prek]][link-prek] |
|   **Meta**    |         [![Hatch project][badge-hatch]][link-hatch] [![Ruff][badge-ruff]][link-ruff] [![uv][badge-uv]][link-uv] [![License][badge-license]][link-license]          |
|  **Package**  |                                                                 [![PyPI][badge-pypi]][link-pypi] [![PyPI][badge-python-versions]][link-pypi]                                                                 |
| **Ecosystem** |                                                                                  [![scverse][badge-scverse]][link-scverse]                                                                                   |
|               |                                                                                                                                                                                                              |

</div>

[badge-scverse]: https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/srivarra/3d16b1f0631693cd529aeabeb089985c/raw/scverse-badge.json
[badge-build]: https://github.com/srivarra/annsel/actions/workflows/build.yaml/badge.svg
[badge-test]: https://github.com/srivarra/annsel/actions/workflows/test.yaml/badge.svg
[badge-docs]: https://img.shields.io/readthedocs/annsel?logo=readthedocs
[badge-codecov]: https://codecov.io/gh/srivarra/annsel/graph/badge.svg?token=ST0ST1BTWU
[badge-ruff]: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
[badge-uv]: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json
[badge-license]: https://img.shields.io/badge/License-MIT-yellow.svg
[badge-hatch]: https://img.shields.io/badge/%F0%9F%A5%9A-Hatch-4051b5.svg
[badge-pypi]: https://img.shields.io/pypi/v/annsel.svg?logo=pypi&label=PyPI&logoColor=gold
[badge-python-versions]: https://img.shields.io/pypi/pyversions/annsel.svg?logo=python&label=Python&logoColor=gold
[badge-prek]: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/j178/prek/master/docs/assets/badge-v0.json
Annsel is a user-friendly library that brings familiar dataframe-style operations to [`AnnData`](https://anndata.readthedocs.io/en/latest/) objects.

It's built on the [narwhals][link-narwhals] compatibility layer for dataframes.

Take a look at the GitHub Projects board for features and future plans: [Annsel Features][link-gh-project]

<!-- done -->

## Getting started

Please refer to the [documentation][link-docs], in particular, the [API documentation][link-api].

There's also a brief tutorial on how to use all the features of `annsel`: [All of Annsel][link-tutorial].

## Installation

You need to have Python 3.12 or newer installed on your system. If you don't have
Python installed, we recommend installing [uv][link-uv].
There are several ways to install `annsel`:

1. Install the most recent release:

    With `uv`:

    ```zsh
    uv add annsel
    ```

    With `pip`:

    ```zsh
    pip install annsel
    ```

2. Install the latest development version:

    With `uv`:

    ```zsh
    uv add git+https://github.com/srivarra/annsel
    ```

    With `pip`:

    ```zsh
    pip install git+https://github.com/srivarra/annsel.git@main
    ```

## Examples

`annsel` comes with a small dataset from Cell X Gene to help you get familiar with the API.

```python
import annsel as an

adata = an.datasets.leukemic_bone_marrow_dataset()
```

The dataset looks like this:

```shell
AnnData object with n_obs × n_vars = 31586 × 458
    obs: 'Cluster_ID', 'donor_id', 'Sample_Tag', 'Cell_label', 'is_primary_data', 'organism_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'assay_ontology_term_id', 'tissue_ontology_term_id', 'Genotype', 'development_stage_ontology_term_id', 'sex_ontology_term_id', 'disease_ontology_term_id', 'cell_type_ontology_term_id', 'suspension_type', 'tissue_type', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid'
    var: 'vst.mean', 'vst.variance', 'vst.variance.expected', 'vst.variance.standardized', 'vst.variable', 'feature_is_filtered', 'Unnamed: 0', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length', 'feature_type'
    uns: 'cell_type_ontology_term_id_colors', 'citation', 'default_embedding', 'schema_reference', 'schema_version', 'title'
    obsm: 'X_bothumap', 'X_pca', 'X_projected', 'X_projectedmean', 'X_tsneni', 'X_umapni'

```

### Filter

You can filter on `obs`, `var`, `var_names`, `obs_names`, `X` and it's layers, as well as `obsm` and `varm` matrices as a key-value pair containing the attribute's key name and the predicate to filter on. *Currently the column names are numerical indices for `obsm` and `varm` matrices.*

```python
adata.an.filter(
    obs=(
        an.col(["Cell_label"]).is_in(["Classical Monocytes", "CD8+CD103+ tissue resident memory T cells"]),
        an.col(["sex"]) == "male",
    ),
    var=an.col(["vst.mean"]) >= 3,
    obsm={"X_pca": an.col([0]) > 0}, # PC1 values greater than 0
    copy=False, # Whether to return a copy of the AnnData object or just a view of it.
)
```

```shell
View of AnnData object with n_obs × n_vars = 736 × 67
    obs: 'Cluster_ID', 'donor_id', 'Sample_Tag', 'Cell_label', 'is_primary_data', 'organism_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'assay_ontology_term_id', 'tissue_ontology_term_id', 'Genotype', 'development_stage_ontology_term_id', 'sex_ontology_term_id', 'disease_ontology_term_id', 'cell_type_ontology_term_id', 'suspension_type', 'tissue_type', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid'
    var: 'vst.mean', 'vst.variance', 'vst.variance.expected', 'vst.variance.standardized', 'vst.variable', 'feature_is_filtered', 'Unnamed: 0', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length', 'feature_type'
    uns: 'cell_type_ontology_term_id_colors', 'citation', 'default_embedding', 'schema_reference', 'schema_version', 'title'
    obsm: 'X_bothumap', 'X_pca', 'X_projected', 'X_projectedmean', 'X_tsneni', 'X_umapni'
```

### Select

You can select on `obs`, `var`, `var_names`, `obs_names`, `X` and it's layers. Selecting returns a new AnnData object. It's useful if you don't need all the columns in `obs` or `var` and just want to work with a few.

```python
adata.an.select(
    obs=an.col(["Cell_label"]),
    var=an.col(["vst.mean", "vst.std"]),
)
```

### Group By

You can group over `obs` and `var` columns which returns a generator of objects containing the grouped data and the grouping parameters.

```python
gb_adata_result = adata.an.group_by(
    obs=an.col(["Cell_label"]),
    var=an.col(["feature_type"]),
    copy=False,
)
```

Here's what the first group looks like:

```python
next(adata.an.group_by(
    obs=an.col(["Cell_label"]),
    copy=False,
))
```

```shell
GroupByAnnData:
  ├── Observations:
  │   └── Cell_label: Lymphomyeloid prog
  ├── Variables:
  │   └── (all variables)
  └── AnnData:
      View of AnnData object with n_obs × n_vars = 913 × 458
          obs: 'Cluster_ID', 'donor_id', 'Sample_Tag', 'Cell_label', 'is_primary_data', 'organism_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'assay_ontology_term_id', 'tissue_ontology_term_id', 'Genotype', 'development_stage_ontology_term_id', 'sex_ontology_term_id', 'disease_ontology_term_id', 'cell_type_ontology_term_id', 'suspension_type', 'tissue_type', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid'
          var: 'vst.mean', 'vst.variance', 'vst.variance.expected', 'vst.variance.standardized', 'vst.variable', 'feature_is_filtered', 'Unnamed: 0', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length', 'feature_type'
          uns: 'cell_type_ontology_term_id_colors', 'citation', 'default_embedding', 'schema_reference', 'schema_version', 'title'
          obsm: 'X_bothumap', 'X_pca', 'X_projected', 'X_projectedmean', 'X_tsneni', 'X_umapni'
```

### Pipe

There's also a small utility method which allows you to chain operations together like in `Xarray` and `Pandas` called `pipe`.

```python
import scanpy as sc
adata.an.pipe(sc.pl.embedding, basis="X_tsneni", color="Cell_label")
```

## Release notes

See the [changelog][changelog].

## Contact

For questions and help requests, you can reach out in the [scverse discourse][scverse-discourse] or the [discussions][link-disucssions] tab.
If you found a bug, please use the [issue tracker][issue-tracker].

## Citation

> Varra, S. R. annsel [Computer software]. <https://github.com/srivarra/annsel>

<!-- done3 -->

[scverse-discourse]: https://discourse.scverse.org/
[issue-tracker]: https://github.com/srivarra/annsel/issues
[changelog]: https://annsel.readthedocs.io/en/latest/changelog.html
[link-docs]: https://annsel.readthedocs.io
[link-api]: https://annsel.readthedocs.io/en/latest/api/index.html
[link-tutorial]: https://annsel.readthedocs.io/en/latest/notebooks/all_of_annsel.html
[link-pypi]: https://pypi.org/project/annsel
[link-codecov]: https://codecov.io/gh/srivarra/annsel
[link-test]: https://github.com/srivarra/annsel/actions/workflows/test.yml
[link-build]: https://github.com/srivarra/annsel/actions/workflows/build.yaml
[link-ruff]: https://github.com/astral-sh/ruff
[link-uv]: https://github.com/astral-sh/uv
[link-license]: https://opensource.org/licenses/MIT
[link-hatch]: https://github.com/pypa/hatch
[link-narwhals]: https://github.com/narwhals-dev/narwhals
[link-disucssions]: https://github.com/srivarra/annsel/discussions
[link-prek]: https://github.com/j178/prek
[link-gh-project]: https://github.com/users/srivarra/projects/9
[link-scverse]: https://scverse.org/packages/#ecosystem
