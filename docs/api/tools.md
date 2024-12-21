# Tools: `tl`

```{eval-rst}
.. module:: annsel.tl
```

```{eval-rst}
.. currentmodule:: annsel
```

## AnnData Accessor: `an`

In order to use the `AnnData` `an` accessor, you need to import `annsel`.

Something like the following:

```python
import annsel
```

or

```python
import annsel as an
```

This will register the `an` accessor to the `AnnData` object.

```{eval-rst}
.. autosummary::
    :toctree: ../generated
    :template: autosummary/accessor_method.rst

    AnnselAccessor.filter
    AnnselAccessor.select
    AnnselAccessor.group_by
    AnnselAccessor.pipe
```
