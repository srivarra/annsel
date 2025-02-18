# annsel

<div align="center">

|               |                                                                                                                                                                                                              |
| :-----------: | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
|  **Status**   | [![Build][badge-build]][link-build] [![Tests][badge-test]][link-test] [![Documentation][badge-docs]][link-docs] [![codecov][badge-codecov]][link-codecov] [![pre-commit][badge-pre-commit]][link-pre-commit] |
|   **Meta**    |         [![Hatch project][badge-hatch]][link-hatch] [![Ruff][badge-ruff]][link-ruff] [![uv][badge-uv]][link-uv] [![License][badge-license]][link-license] [![gitmoji][badge-gitmoji]][link-gitmoji]          |
|  **Package**  |                                                                 [![PyPI][badge-pypi]][link-pypi] [![PyPI][badge-python-versions]][link-pypi]                                                                 |
| **Ecosystem** |                                                                                  [![scverse][badge-scverse]][link-scverse]                                                                                   |
|               |                                                                                                                                                                                                              |

</div>

[badge-scverse]: https://camo.githubusercontent.com/67057d1bc0b82681113a209b6478a92476a70ed2876c25de957a2090e04cb0bd/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f736376657273652d65636f73797374656d2d626c75652e7376673f6c6162656c436f6c6f723d7768697465266c6f676f3d646174613a696d6167652f737667253262786d6c3b6261736536342c5044393462577767646d567963326c76626a30694d5334774969426c626d4e765a476c755a7a3069565652474c54676949484e305957356b59577876626d5539496d3576496a382b5043464554304e555756424649484e325a79425156554a4d53554d67496930764c31637a517938765246524549464e57527941784c6a45764c30564f496941696148523063446f764c336433647935334d793576636d637652334a686347687059334d7655315a484c7a45754d5339455645517663335a6e4d5445755a48526b496a344b50484e325a79423361575230614430694d5441774a534967614756705a326830505349784d44416c4969423261575633516d393450534977494441674f5445674f54456949485a6c636e4e7062323439496a45754d53494b49434167494868746247357a50534a6f644852774f693876643364334c6e637a4c6d39795a7938794d4441774c334e325a79494b49434167494868746247357a4f6e68736157357250534a6f644852774f693876643364334c6e637a4c6d39795a7938784f546b354c3368736157357249694234625777366333426859325539496e42795a584e6c636e5a6c49676f674943416765473173626e4d366332567961575939496d6830644841364c793933643363756332567961575975593239744c7949676333523562475539496d5a7062477774636e56735a54706c646d56756232526b4f324e7361584174636e56735a54706c646d56756232526b4f334e30636d39725a5331736157356c616d3970626a7079623356755a44747a64484a766132557462576c305a584a736157317064446f794f79492b43694167494341385a7942705a44306952574a6c626d56664d79492b4369416749434167494341675047632b43694167494341674943416749434167494478775958526f49475139496b307a4e5377344f533432597930794d69347a4c43307a4c6a51674c544d774c6a59734c5445354c6a67674c544d774c6a59734c5445354c6a686a4d5441754f4377784e6934354944517a4c446b754d5341314d6934354c4449754e574d784d6934304c4330344c6a4d674f4377744d5455754d7941324c6a67734c5445344c6a466a4e5334304c4463754d6941314c6a4d734d6a4d754e5341744d5334784c4449354c6a526a4c5455754e6977314c6a45674c5445314c6a4d734e793435494330794f437732576949676333523562475539496d5a706247773649325a6d5a6a746d615778734c584a3162475536626d3975656d5679627a747a64484a7661325536497a41774d44747a64484a766132557464326c6b644767364d5842344f79497650676f6749434167494341674943416749434138634746306143426b50534a4e4f444d754f5377304d793431597a49754f5377744e793478494441754f4377744d5449754e5341774c6a55734c54457a4c6a4e6a4c5441754e7977744d53347a494330784c6a55734c5449754d7941744d6934304c43307a4c6a466a4c5445324c6a45734c5445794c6a59674c5455314c6a6b734d5341744e7a41754f5377784e693434597930784d4334354c4445784c6a55674c5445774c6a45734d6a41674c5459754e7977794e533434597a4d754d5377304c6a67674e7934354c4463754e6941784d7934304c446c6a4c5445784c6a55734c5445794c6a51674f5334344c43307a4d533478494449354c43307a4f474d794d5377744e79343149444d794c6a55734c544d674d7a63754d5377794c6a68614969427a64486c735a5430695a6d6c7362446f6a4d7a517a4e444d304f325a7062477774636e56735a547075623235365a584a764f334e30636d39725a546f6a4d4441774f334e30636d39725a5331336157523061446f78634867374969382b43694167494341674943416749434167494478775958526f49475139496b30334f5334324c4455774c6a526a4f5377744d5441754e5341314c4330784f533433494451754f4377744d6a41754e474d744d437777494451754e4377334c6a45674d6934794c4449794c6a5a6a4c5445754d6977344c6a55674c5455754e4377784e6941744d5441754d5377784d533434597930794c6a45734c5445754f4341744d7977744e69343549444d754d5377744d5452614969427a64486c735a5430695a6d6c7362446f6a5a6d5a6d4f325a7062477774636e56735a547075623235365a584a764f334e30636d39725a546f6a4d4441774f334e30636d39725a5331336157523061446f78634867374969382b43694167494341674943416749434167494478775958526f49475139496b30324e4377314e4334795979307a4c6a4d734c5451754f4341744f4334784c4330334c6a51674c5445794c6a4d734c5445774c6a686a4c5449754d6977744d533433494330784e6934304c4330784d533479494330784f5334794c4330784e533478597930324c6a51734c5459754e4341744f5334314c4330784e6934354943307a4c6a51734c54497a4c6a466a4c5451754e4377744d433434494330344c6a49734d433479494330784d4334324c4445754e574d744d5334784c4441754e6941744d6934784c4445754d6941744d6934344c444a6a4c5459754e7977324c6a49674c5455754f4377784e7941744d5334324c4449304c6a4e6a4e4334314c4463754f4341784d7934794c4445314c6a51674d6a51754d7977794d693434597a55754d53777a4c6a51674d5455754e6977344c6a51674d546b754d7977784e6d4d784d5334334c4330344c6a45674e7934324c4330784e433435494459754d7977744d5463754e6c6f6949484e306557786c50534a6d615778734f694e694e474930596a51375a6d6c73624331796457786c4f6d3576626e706c636d3837633352796232746c4f694d774d444137633352796232746c4c5864705a48526f4f6a4677654473694c7a344b4943416749434167494341674943416750484268644767675a44306954544d344c6a63734f533434597a63754f5377324c6a4d674d5449754e4377354c6a67674d6a41734f433431597a55754e7977744d5341304c6a6b734c5463754f5341744e4377744d544d754e6d4d744e4334304c4330794c6a67674c546b754e4377744e433479494330784e5334334c4330304c6a4a6a4c5463754e5377744d4341744d5459754d79777a4c6a6b674c5449774c6a59734e693430597a51734c5449754d7941784d5334354c43307a4c6a67674d6a41754d7977794c6a6c614969427a64486c735a5430695a6d6c7362446f6a5a6d5a6d4f325a7062477774636e56735a547075623235365a584a764f334e30636d39725a546f6a4d4441774f334e30636d39725a5331336157523061446f78634867374969382b4369416749434167494341675043396e50676f67494341675043396e50676f384c334e325a7a343d
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
[badge-pre-commit]: https://results.pre-commit.ci/badge/github/srivarra/annsel/main.svg
[badge-gitmoji]: https://img.shields.io/badge/gitmoji-😜😍-FFDD67.svg

Annsel is a user-friendly library that brings familiar dataframe-style operations to [`AnnData`](https://anndata.readthedocs.io/en/latest/) objects.

It's built on the [narwhals][link-narwhals] compatibility layer for dataframes.

Take a look at the GitHub Projects board for features and future plans: [Annsel Features][link-gh-project]

<!-- done -->

## Getting started

Please refer to the [documentation][link-docs], in particular, the [API documentation][link-api].

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

```python
import annsel as an

adata=an.datasets.leukemic_bone_marrow_dataset()
```

### Filter

```python
adata.an.filter(
    obs=(
        an.col(["Cell_label"]).is_in(["Classical Monocytes", "CD8+CD103+ tissue resident memory T cells"]),
        an.col(["sex"]) == "male",
    ),
    var=an.col(["vst.mean"]) >= 3,
)
```

### Select

```python
adata.an.select(
    obs=an.col(["Cell_label"]),
    var=an.col(["vst.mean", "vst.std"]),
)
```

### Group By

```python
adata.an.group_by(
    obs=an.col(["Cell_label"]),
    var=an.col(["feature_type"]),
    return_group_names=True,
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
[link-pre-commit]: https://results.pre-commit.ci/latest/github/srivarra/annsel/main
[link-gitmoji]: https://gitmoji.dev/
[link-gh-project]: https://github.com/users/srivarra/projects/9
[link-scverse]: https://scverse.org/packages/#ecosystem
