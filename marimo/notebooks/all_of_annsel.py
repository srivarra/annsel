# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "narwhals==1.22",
#     "marimo",
#     "anndata==0.11.3",
#     "annsel==0.0.4",
# ]
# ///
import marimo

__generated_with = "0.10.13"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# All of Annsel""")
    return


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        This notebook tells you just about everything you need to use `annsel`. It's a good starting point to get a feel for the package.

        :::{note}
        :class: dropdown

        You should be familiar with [`AnnData` ](https://anndata.readthedocs.io/en/latest/) beforehand.
        :::
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Set up Data""")
    return


@app.cell
def _():
    # '%load_ext autoreload\n%autoreload 2' command supported automatically in marimo
    return


@app.cell
def _():
    import annsel as an

    return (an,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""We will load a Leukemic bone marrow cytometry dataset :cite:p:`Triana2021` You can view the dataset on [CellxGene](https://cellxgene.cziscience.com/e/b3a5a10f-b1cb-4e8e-abce-bf345448625b.cxg/)."""
    )
    return


@app.cell
def _(an):
    adata = an.datasets.leukemic_bone_marrow_dataset()
    return (adata,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""View the contents of the `AnnData` object.""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        Importing `annsel` will automatically register the accessor and add the `an` attribute to the `AnnData` object's namespace.

        You can access the methods of the accessor using the `an` attribute on an `AnnData` object.

        ```python

        anndata.AnnData.an.<method_name>
        ```
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Filter


        We can filter the `AnnData` object using the `filter` method. Filtering can be applied to the `obs`, `var`, `X` (with a layer), `obs_names` and `var_names`.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Let's filter the `AnnData` object by the Cell Type stored in `obs["Cell_label"]`.""")
    return


@app.cell
def _(adata, an):
    adata.an.filter(obs=an.col(["Cell_label"]) == "CD8+CD103+ tissue resident memory T cells")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        :::{note}
        This is equivalent to

        ```python
        adata[adata.obs["Cell_label"] == "CD8+CD103+ tissue resident memory T cells", :]
        ```
        :::
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""You can also apply other expressions to `an.col` callable, for example `is_in` allows you to filter by a list of values. Under the hood, `annsel` uses [Narwhals](https://narwhals-dev.github.io/narwhals/) to apply these expressions. A full list of expressions which can be applied to a column can be found [here](https://narwhals-dev.github.io/narwhals/api-reference/expr/)."""
    )
    return


@app.cell
def _(adata, an):
    adata.an.filter(
        obs=an.col(["Cell_label"]).is_in(["Classical Monocytes", "CD8+CD103+ tissue resident memory T cells"])
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""We can also combine multiple Predicates using the `&` and `|` operators.""")
    return


@app.cell
def _(adata, an):
    adata.an.filter(
        obs=(an.col(["Cell_label"]).is_in(["Classical Monocytes", "CD8+CD103+ tissue resident memory T cells"]))
        & (an.col(["sex"]) == "male")
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Or if you pass a tuple of Predicates, it will apply the `&` operator between them automatically.""")
    return


@app.cell
def _(adata, an):
    adata.an.filter(
        obs=(
            an.col(["Cell_label"]).is_in(["Classical Monocytes", "CD8+CD103+ tissue resident memory T cells"]),
            an.col(["sex"]) == "male",
        )
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""The `|` operator is applied between the Predicates as well. Here we will select all the cells that are either Classical Monocytes or CD8+CD103+ tissue resident memory T cells, or are from male samples irrespective of the cell type."""
    )
    return


@app.cell
def _(adata, an):
    adata.an.filter(
        obs=(
            an.col(["Cell_label"]).is_in(["Classical Monocytes", "CD8+CD103+ tissue resident memory T cells"])
            | (an.col(["sex"]) == "male")
        )
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""We can also filter the `AnnData` object by the `var` column. Here we will filter the `AnnData` object to only include the genes with `vst.mean` greater than or equal to 3."""
    )
    return


@app.cell
def _(adata, an):
    adata.an.filter(var=an.col(["vst.mean"]) >= 3)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""Filtering can also be applied to the `X` matrix. Here we will filter the `AnnData` object to only include the cells with `ENSG00000205336` gene expression greater than 1. If you want to filter by a layer, you can pass the layer name to the `layer` argument. and the operation will be applied to the `var_name` of that layer."""
    )
    return


@app.cell
def _(adata, an):
    adata.an.filter(
        x=an.col(["ENSG00000205336"]) > 1,
        layer=None,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""This can all be combined together as well""")
    return


@app.cell
def _(adata, an):
    adata.an.filter(
        obs=(
            an.col(["Cell_label"]).is_in(["Classical Monocytes", "CD8+CD103+ tissue resident memory T cells"]),
            an.col(["sex"]) == "male",
        ),
        var=an.col(["vst.mean"]) >= 3,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""Filtering can also be applied to `var_names` and `obs_names` using the `an.var_names` and `an.obs_names` predicates respectively. These are special predicates which can only be applied to `var_names` and `obs_names` of the `AnnData` object."""
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""Here, arbitrarily, we will filter the AnnData object to only include the cells with `obs_names` starting with `645` and the genes with `var_names` starting with `ENSG0000018`."""
    )
    return


@app.cell
def _(adata, an):
    adata.an.filter(
        obs_names=an.obs_names.str.starts_with("645"), var_names=an.var_names.str.starts_with("ENSG0000018")
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Select


        We can also apply the `select` method to the `AnnData` object. This is similar to the `filter` method, but it will only keep the rows and columns that match the Predicates. It can be applied to the `obs`, `var`, and `X`.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""Here we will select the `Cell_label` and `sex` columns from the `obs` table, the `feature_name` column from the `var` table and the `ENSG00000205336` gene from `X`. This will return a new `AnnData` object with only these columns in the `obs`, `var` and `X` tables."""
    )
    return


@app.cell
def _(adata, an):
    adata.an.select(
        obs=an.col(["Cell_label", "sex"]),
        var=an.col(["feature_name"]),
        x=an.col(["ENSG00000205336"]),
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Group By


        We can also group the `AnnData` object by the `obs` and `var` columns. This will return a generator of `AnnData` objects subset on each group.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""Here we will group the `AnnData` object by the `Cell_label` column in the `obs` table and the `feature_type` column in the `var` table."""
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""If you pass `return_group_names=True`, the generator will yield a tuple of the group name and the `AnnData` object. If you group by both `obs` and `var`, the generator will yield a tuple of the both group names and the `AnnData` object, if you group by only one, it will yield a tuple of the group name and the `AnnData` object."""
    )
    return


@app.cell
def _(adata, an):
    for i in adata.an.group_by(
        obs=an.col(["Cell_label"]),
        var=an.col(["feature_type"]),
        return_group_names=True,
    ):
        obs_group, var_group, _adata = i
    return i, obs_group, var_group


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
