# API and mathematical conventions

## TDatum

```python
from tdatum import TDatum
td = TDatum(A_plus, A_minus, D='identity', check=True)
```

The inputs are nonempty square Sage matrices over a univariate Laurent
polynomial ring. They must use the same variable and have integer coefficients.
The T-datum conditions force polynomial support in nonnegative degrees.
`D` is a positive integer diagonal matrix of the same size.

The constructor recovers `N0`, `N_plus`, and `N_minus` by coefficientwise
positive and negative parts. It checks the common positive part,
the permutation-monomial form of `N0`, positive degrees, nonnegative
coefficients, interior support, and disjoint support. These are the
conditions (N1)--(N4) in the
[mathematical reference](https://arxiv.org/abs/1912.05710).

By default it also checks `N0 D = D N0`, integrality of
`D^(-1) N_plus D` and `D^(-1) N_minus D`, and the symplectic relation
`A_plus D A_minus^dagger = A_minus D A_plus^dagger`.
Here `A^dagger` means substitution of `z^(-1)` for `z` followed by transpose.
`check=False` skips only these three compatibility checks and is intended
for use when the identities are already known.

| Method | Result |
| --- | --- |
| `pair()` | `(A_plus, A_minus)` |
| `triple()` | `(N0, N_plus, N_minus)` |
| `degrees()` | Tuple of positive row degrees |
| `symmetrizer()` | Diagonal of `D` as a tuple |
| `permutation()` | Sage permutation in the leading monomials of `N0` |
| `sign_dual()` | `(A_minus, A_plus, D)` as a T-datum |
| `langlands_dual()` | Conjugated matrices and dual symmetrizer |
| `mutation_loop()` | Mutation loop on the maximal initial index set |
| `maximal_initial_indices()` | Zero-based pairs `(a,p)` with `0 <= p < p_a` |
| `is_indecomposable()` | Whether the datum is indecomposable |
| `connected_component(a)` | Initial indices in the component through `(a,0)` |
| `plot_mutation_loop()` | Sage graphics for the maximal initial quiver |

`mutation_loop(R=...)` accepts a consistent subset of the initial indices,
in the sense of the reference. The caller is responsible for choosing such
a subset; the resulting mutation loop is checked by its constructor.
The stored matrices are immutable. Use `matrix(td.pair()[0])` to obtain a
mutable copy.

### Plot a quiver

```python
from tdatum.examples import RSG

td = RSG([3, 1]).t_datum()
td.plot_mutation_loop().show(figsize=6, axes_pad=0.15)
```

The plot labels vertices by the zero-based pairs `(a, p)` in
`maximal_initial_indices()` and highlights the first mutation block (`p = 0`)
in green. To plot the exchange matrix at time `u` with SageMath's standard
quiver layout, use:

```python
from sage.all import ClusterQuiver

loop = td.mutation_loop()
u = 1
ClusterQuiver(loop.b_matrix(u)).plot().show(figsize=6, axes_pad=0.15)
```

Here the vertex labels are the integer indices of the exchange matrix.
See the [RSG notebook](../examples/02_rsg.ipynb) for the pair-labeled plot
and the [mutation-loop notebook](../examples/03_mutation_loops.ipynb) for
plots at successive time steps with one-based display labels.

## MutationLoop

```python
from tdatum import MutationLoop
loop = MutationLoop(B, sequence_of_indices, permutation, symmetrizer=None)
```

`B` is an integer exchange matrix and `symmetrizer` is a tuple of positive
integers `d` such that `B * diagonal_matrix(d)` is skew-symmetric.
The default is the tuple of ones. Each block in `sequence_of_indices`
contains distinct zero-based vertices that commute at that time.
The permutation uses Sage's one-based convention, preserves `d`, and
identifies the endpoint matrix with the initial matrix.

| Method | Result |
| --- | --- |
| `vertices()` | Copy of the mutation blocks in one loop |
| `indices(u)` | Copy of the block at a time step within one loop |
| `length()` | Number of time steps in one loop |
| `whole_length()` | Total number of mutations in one loop |
| `size()` | Number of quiver vertices |
| `initial_b_matrix()` | Initial exchange matrix |
| `b_matrix(u)` | Exchange matrix at nonnegative integer time `u` |
| `permutation()` | Endpoint permutation |
| `symmetrizer()` | Right symmetrizer |
| `inverse()` | Inverse loop with the same symmetrizer |
| `is_complete()` | Whether every vertex orbit meets a mutation point |
| `t_datum(variable_name='z')` | Validated T-datum of a complete loop |

Equality compares the initial matrix, symmetrizer, endpoint permutation,
and time blocks. The order of commuting vertices inside one block is
ignored. Time blocks themselves are retained.

The three notions of time length, number of mutations, and nonlinear period
are distinct. The library does not determine a nonlinear period.

## Example constructors

Import constructors from `tdatum.examples`. Every constructor has `pair()`
and `t_datum(D='identity')`. `pair()` constructs matrices;
`t_datum()` additionally runs the validation described above.

| Constructor | Scope of the examples tested for this development version |
| --- | --- |
| `RSG(n_list)` | `[2,1]`, `[3,1]`, `[4,1]`, `[3,2,2]`, `[4,2,3]`, `[3,2,2,2]` |
| `Rank2(label)` | All six stored labels `1,...,6` |
| `LengthOne(degree_list)` | Coefficient list `[1,-2,1]` |
| `SG(n_list)` | `[3,2]` |
| `UntwistedAffine(type, rank, level)` | `('A',2,2)` |
| `TamelyLaced(...)`, `Unknown(...)` | Legacy constructors retained; no new release verification |

For `RSG`, the first parameter is at least 2, all parameters are positive
integers, and their sum exceeds 2. The input-domain check does not establish
validity for every list; use `t_datum()` to validate a concrete construction.
The table records finite tests, not a classification or an all-parameter
proof. The six rank-two labels are a stored catalogue; the package makes
no exhaustiveness claim.
