# TDatum for SageMath

Compute with T-data and their mutation loops in SageMath.

`TDatum` checks a pair of polynomial matrices and a symmetrizer, recovers
their defining data, constructs mutation loops, and computes sign and
Langlands duals. `MutationLoop` constructs T-data from complete mutation
loops. The example catalogue includes reduced sine-Gordon and rank-two
examples.

The mathematical reference is Yuma Mizuno,
[Difference equations arising from cluster algebras](https://arxiv.org/abs/1912.05710).

This is development version `0.1.0.dev0`, tested with SageMath 10.8.
The [GitHub repository](https://github.com/yuma-mizuno/tdatum) is currently
private while release preparation is in progress.

## Install

First install SageMath using its [installation guide](https://doc.sagemath.org/html/en/installation/).
Clone the repository using a GitHub account with access, then install into
your SageMath environment:

```sh
git clone https://github.com/yuma-mizuno/tdatum.git
cd tdatum
sage -pip install --no-deps .
```

For development use `sage -pip install --no-deps -e .`.
SageMath is a prerequisite: installing this package does not install SageMath.
The code examples below run in a SageMath notebook or with `sage -python`.
They also use ordinary Python exponentiation (`**`), so their meaning does not
depend on the Sage preparser.

## First example

```python
from sage.all import QQ, LaurentPolynomialRing, matrix
from tdatum import TDatum

R = LaurentPolynomialRing(QQ, "z")
z = R.gen()
A_plus = matrix(R, [[1 + z**2, 0], [0, 1 + z**2]])
A_minus = matrix(R, [[1 + z**2, -z], [-z, 1 + z**2]])
td = TDatum(A_plus, A_minus)

assert td.degrees() == (2, 2)
assert td.symmetrizer() == (1, 1)
loop = td.mutation_loop()
assert loop.is_complete()
assert loop.t_datum().pair() == td.pair()
```

Construct a reduced sine-Gordon example:

```python
from tdatum.examples import RSG

td = RSG([3, 1]).t_datum()
assert td.degrees() == (2, 6)
```

Plot its initial quiver:

```python
td.plot_mutation_loop().show(figsize=6, axes_pad=0.15)
```

Vertices are labeled by the zero-based pairs `(a, p)` in
`td.maximal_initial_indices()`. Green vertices have `p = 0` and form the
first mutation block. The [mutation-loop notebook](examples/03_mutation_loops.ipynb)
also plots the quivers before and after each time step.

`RSG(...).pair()` returns the constructed matrices. `RSG(...).t_datum()`
also validates them. The same interface is available on the other example
constructors; pass `D` explicitly when a nonidentity symmetrizer is needed.

## Examples and conventions

- [Getting started](examples/01_getting_started.ipynb): construction,
  nonidentity symmetrizers, and duals.
- [RSG examples](examples/02_rsg.ipynb): six explicit parameter lists,
  their matrices, degrees, and an initial quiver plot.
- [Mutation loops](examples/03_mutation_loops.ipynb): quivers along a loop,
  time steps, vertex labels, endpoint permutations, and reconstruction of a T-datum.
- [API and mathematical conventions](docs/api.md).
- [Verification record](docs/verification.md).

Python matrix indices and mutation vertices are zero-based. Sage permutations
are one-based. The notebooks display mathematical vertex labels starting at one
and translate to the API at the construction step.

The default constructor checks the T-datum conditions, including integrality,
(N1)--(N4), and the compatibility identities for `D`. This is an exact check of
the supplied finite input. A T-datum need not have a periodic T/Y-system;
finite type and periodicity require separate mathematical arguments.

## Test

```sh
sage -python tools/check.py
```

This runs regression tests and Sage doctests. To rebuild and execute the
notebooks and export their HTML versions:

```sh
sage -python tools/build_examples.py
```

Notebook generation uses `nbformat`, `nbclient`, and `nbconvert`; the runtime
library itself uses SageMath and Python's standard library.

## Development and citation

See [CONTRIBUTING.md](CONTRIBUTING.md) for contributions and
[CITATION.cff](CITATION.cff) for citation metadata. The code is distributed
under GPL-3.0-or-later; see [LICENSE](LICENSE). Source history and the public
release preparation are recorded in [docs/provenance.md](docs/provenance.md) and
[docs/release-preparation.md](docs/release-preparation.md).
