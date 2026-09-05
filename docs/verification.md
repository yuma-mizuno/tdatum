# Verification record

Verified locally on 2026-09-05 with SageMath 10.8 and its Python 3.12
environment in WSL.

| Check | Result |
| --- | --- |
| Regression suite | 12 test methods passed, including parameterized invalid inputs and example families |
| Sage doctests | 301 examples passed, fixed random seed 0 |
| Wheel build | `sagemath_tdatum-0.1.0.dev0-py3-none-any.whl` built |
| Isolated installation | Wheel imported from a temporary installation; all 12 regression tests passed there |
| Getting-started notebook | 5 code cells executed, 0 errors |
| RSG notebook | 5 code cells executed, 0 errors; initial quiver plot embedded |
| Mutation-loop notebook | 6 code cells executed, 0 errors; quivers at three time steps embedded |
| Citation metadata | Validated against the official CFF 1.2.0 JSON schema |

## GitHub Actions

The first full [GitHub Actions run](https://github.com/yuma-mizuno/tdatum/actions/runs/33951789540)
passed on 2026-09-05 for commit `2148ecec96abec857e3b55938f670ff7ec298c07`.
The SageMath 10.8 container passed package installation, all 12 regression
test methods, all 301 doctests, the isolated wheel installation, and all
14 code cells across the three notebooks. Validation logs and notebook HTML
exports are attached to that run as an artifact.

Subsequent results are available in the
[workflow history](https://github.com/yuma-mizuno/tdatum/actions/workflows/tests.yml).
Repository access is required while the project is private.

## Reproduce

From the repository root:

```sh
sage -python tools/check.py
sage -python tools/check_install.py
sage -python tools/build_examples.py
```

The test tools write logs to `test-results/`. Notebook execution updates the
three files in `examples/` and exports HTML to `docs/_build/examples/`.
Those logs and HTML exports are ignored by Git; the executed notebooks are
tracked. GitHub Actions runs the same commands and uploads validation
artifacts with a retention period of 14 days.

## Properties covered

The regression suite covers malformed symmetrizers, noninteger coefficients,
inconsistent positive parts, nonpositive degrees, overlapping supports,
support outside the permitted degree interval, and compatibility of `D`.
It checks sign and Langlands duality, preservation of a nonidentity
symmetrizer through loop construction and inversion, a nonsymmetric
datum-to-loop-to-datum round trip, and custom polynomial variable names.

For mutation loops, it checks the distinction between different time blocks,
invariance under reordering within a commuting block, validation of blocks
and symmetrizers, conversion of complete loops, and protection from mutation
of input containers.

The representative families are listed in [api.md](api.md). The RSG matrix
pairs for the six listed inputs were also checked against an independent
coefficient-level implementation of the defining conditions in the local
pre-extraction audit. They were compared exactly with the continued-fraction
implementation for the same inputs.

## Limits of this verification

This record establishes the outcomes of the specified computations. It does
not prove correctness for every possible input, finite type, nonlinear
periodicity, or a classification of T-data.

`TamelyLaced` and `Unknown` are inherited constructors retained for further
review. They are not used by the introductory notebooks and have no new
release verification. Only the listed `SG` and `UntwistedAffine` inputs have
been tested; other parameters and other Sage versions remain to be checked.

The original core had three failing doctests: a nonidentity-symmetrizer
example omitted `D`, and a plot example expected a digraph instead of the
graphics object returned by the implementation. These examples were corrected.
Each standalone doctest now imports the public package explicitly.
