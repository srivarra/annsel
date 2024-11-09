# annsel

<div align="center">

|            |                                                                                                                                                           |
| :--------: | :-------------------------------------------------------------------------------------------------------------------------------------------------------: |
| **Status** | [![Build][badge-build]][link-build] [![Tests][badge-test]][link-test] [![Documentation][badge-docs]][link-docs] [![codecov][badge-codecov]][link-codecov] |
|  **Meta**  |     [![Hatch project][badge-hatch]][link-hatch] [![Ruff][badge-ruff]][link-ruff] [![uv][badge-uv]][link-uv] [![License][badge-license]][link-license]     |
|            |                                                                                                                                                           |

</div>

[badge-build]: https://github.com/srivarra/annsel/actions/workflows/build.yaml/badge.svg
[badge-test]: https://github.com/srivarra/annsel/actions/workflows/test.yaml/badge.svg
[badge-docs]: https://img.shields.io/readthedocs/annsel?logo=readthedocs
[badge-codecov]: https://codecov.io/gh/srivarra/annsel/graph/badge.svg?token=ST0ST1BTWU
[badge-ruff]: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
[badge-uv]: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json
[badge-license]: https://img.shields.io/badge/License-MIT-yellow.svg
[badge-hatch]: https://img.shields.io/badge/%F0%9F%A5%9A-Hatch-4051b5.svg

`annsel` brings familiar DataFrame-style operations to [`AnnData`](https://anndata.readthedocs.io/en/latest/) objects, making filtering and selection intuitive and straightforward. Built on the [narwhals][link-narwhals] library, it provides a seamless interface for manipulating complex biological datasets stored in `AnnData` format.

> [!WARNING]
> This is a very early work in progress, there is set-in-stone API and breaking changes are likely to occur.

## Getting started

Please refer to the [documentation][link-docs],
in particular, the [API documentation][link-api].

## Installation

You need to have Python 3.10 or newer installed on your system. If you don't have
Python installed, we recommend installing [uv][link-uv].
There are several alternative options to install `annsel`:

1. Install the latest development version:

    With `pip`:

    ```zsh
    pip install git+https://github.com/srivarra/annsel.git@main
    ```

    With `uv`:

    ```zsh
    uv add git+https://github.com/srivarra/annsel
    ```

## Release notes

See the [changelog][changelog].

## Contact

For questions and help requests, you can reach out in the [scverse discourse][scverse-discourse].
If you found a bug, please use the [issue tracker][issue-tracker].

## Citation

> t.b.a

[scverse-discourse]: https://discourse.scverse.org/
[issue-tracker]: https://github.com/srivarra/annsel/issues
[changelog]: https://annsel.readthedocs.io/latest/changelog.html
[link-docs]: https://annsel.readthedocs.io
[link-api]: https://annsel.readthedocs.io/en/latest/api.html
[link-pypi]: https://pypi.org/project/annsel
[link-codecov]: https://codecov.io/gh/srivarra/annsel
[link-test]: https://github.com/srivarra/annsel/actions/workflows/test.yml
[link-build]: https://github.com/srivarra/annsel/actions/workflows/build.yaml
[link-ruff]: https://github.com/astral-sh/ruff
[link-uv]: https://github.com/astral-sh/uv
[link-license]: https://opensource.org/licenses/MIT
[link-hatch]: https://github.com/pypa/hatch
[link-narwhals]: https://github.com/narwhals-dev/narwhals
