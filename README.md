# annsel

<div align="center">

<table>
  <tr>
    <th align="right">Package</th>
    <td align="center">
      <a href="https://pypi.org/project/annsel"><img src="https://img.shields.io/pypi/v/annsel?logo=pypi&labelColor=1C2C2E&color=C96329" alt="PyPI"/></a>
      <a href="https://pypi.org/project/annsel"><img src="https://img.shields.io/pypi/pyversions/annsel?logo=python&labelColor=1C2C2E&color=3776AB" alt="Python"/></a>
      <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-yellow?labelColor=1C2C2E" alt="License"/></a>
    </td>
  </tr>
  <tr>
    <th align="right">CI/CD</th>
    <td align="center">
      <a href="https://github.com/srivarra/annsel/actions/workflows/test.yaml"><img src="https://img.shields.io/github/actions/workflow/status/srivarra/annsel/test.yaml?branch=main&logo=github&label=tests&labelColor=1C2C2E" alt="Tests"/></a>
      <a href="https://github.com/srivarra/annsel/actions/workflows/build.yaml"><img src="https://img.shields.io/github/actions/workflow/status/srivarra/annsel/build.yaml?branch=main&logo=github&label=build&labelColor=1C2C2E" alt="Build"/></a>
      <a href="https://codecov.io/gh/srivarra/annsel"><img src="https://codecov.io/gh/srivarra/annsel/graph/badge.svg?token=ST0ST1BTWU" alt="codecov"/></a>
      <a href="https://codspeed.io/srivarra/annsel"><img src="https://img.shields.io/endpoint?url=https://codspeed.io/badge.json" alt="CodSpeed"/></a>
      <a href="https://results.pre-commit.ci/latest/github/srivarra/annsel/main"><img src="https://results.pre-commit.ci/badge/github/srivarra/annsel/main.svg" alt="pre-commit.ci"/></a>
    </td>
  </tr>
  <tr>
    <th align="right">Docs / Other</th>
    <td align="center">
      <a href="https://annsel.readthedocs.io"><img src="https://img.shields.io/readthedocs/annsel?logo=readthedocs&labelColor=1C2C2E" alt="Docs"/></a>
      <a href="https://scverse.org/packages/#ecosystem"><img src="https://img.shields.io/endpoint?url=https%3A%2F%2Fgist.githubusercontent.com%2Fsrivarra%2F3d16b1f0631693cd529aeabeb089985c%2Fraw%2Fscverse-badge.json" alt="scverse"/></a>
    </td>
  </tr>
  <tr>
    <th align="right">Tools</th>
    <td align="center">
      <a href="https://github.com/astral-sh/ruff"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json" alt="Ruff"/></a>
      <a href="https://github.com/astral-sh/uv"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json" alt="uv"/></a>
      <a href="https://gitmoji.dev/"><img src="https://img.shields.io/badge/gitmoji-😜😍-FFDD67?labelColor=1C2C2E" alt="gitmoji"/></a>
    </td>
  </tr>
</table>

</div>

Annsel is a user-friendly library that brings familiar dataframe-style operations to [`AnnData`](https://anndata.readthedocs.io/en/latest/) objects.

It's built on the [narwhals][link-narwhals] compatibility layer for dataframes.

<!-- done -->

## Getting started

Please refer to the [documentation][link-docs], in particular, the [API documentation][link-api].

There's also a brief tutorial on how to use all the features of `annsel`: [All of Annsel][link-tutorial].

## Installation

You need to have Python 3.10 or newer installed on your system. If you don't have
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

```
31,586 cells × 458 genes with cell type annotations and QC metrics
```

### Filter observations and variables

```python
# Filter cells and genes simultaneously
filtered = adata.an.filter(
    obs=an.col("Cell_label").is_in(["Classical Monocytes", "T cells"]),
    var=an.col("vst.mean") >= 3,
)
# Returns: AnnData with filtered cells and genes
```

### Select columns

```python
# Select specific metadata columns
selected = adata.an.select(
    obs=an.col("Cell_label", "sex"),
    var=an.col("feature_name"),
)
# Returns: AnnData with selected columns
```

### Cross-component aggregation

Combine cell metadata with gene expression for powerful analysis:

```python
# Group by cell type, compute marker gene statistics
stats = (
    adata.an.with_obs(
        "Cell_label",
        var_names=["ENSG00000204472", "ENSG00000206560"]
    )
    .group_by("Cell_label")
    .agg([
        nw.col("ENSG00000204472").mean().alias("gene1_mean"),
        nw.col("ENSG00000206560").mean().alias("gene2_mean"),
        nw.len().alias("n_cells"),
    ])
)
# Returns: DataFrame with expression statistics per cell type
```

### Chain operations

```python
# Fully chainable API
result = (
    adata.an
    .filter(obs=an.col("sex") == "male")
    .an.select(var=an.col("feature_name", "vst.mean"))
)
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
[link-uv]: https://github.com/astral-sh/uv
[link-narwhals]: https://github.com/narwhals-dev/narwhals
[link-disucssions]: https://github.com/srivarra/annsel/discussions
