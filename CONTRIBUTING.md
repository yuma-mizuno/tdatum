# Contributing

Use a SageMath environment and install this checkout in editable mode:

```sh
sage -pip install --no-deps -e .
sage -python tools/check.py
```

For a bug report, include the SageMath version, a short input that reproduces
the problem, the expected mathematical property, and the observed result.
Use exact arithmetic where possible.

For a new constructor, document its parameter domain and symmetrizer, cite
the underlying construction, and add representative tests. Keep a matrix-pair
construction separate from the validation of a T-datum and from claims about
periodicity or finite type.

Changes to the mathematical core should include regression tests for the
affected property and run the Sage doctests. Rebuild the notebooks after
changing a documented API. Use UTF-8 for text files.

The maintainer reviews changes on a best-effort basis. A public issue tracker
will be linked here when the repository is published. Contributions are
distributed under the same GPL-3.0-or-later license as the project.
