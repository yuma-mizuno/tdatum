# GitHub release preparation

The project is being prepared in the private repository
[`yuma-mizuno/tdatum`](https://github.com/yuma-mizuno/tdatum).
The distribution name is `sagemath-tdatum`, the Python import name is
`tdatum`, and the development version is `0.1.0.dev0`.

## Scope

The GitHub repository contains the SageMath library, three executed
introductory notebooks, API documentation, tests, build tools, a license,
and citation metadata. English is the primary language for documentation
and examples. GitHub is the distribution and collaboration platform for
this preparation stage.

The introductory sequence is installation, a minimal T-datum, a nonidentity
symmetrizer, RSG examples, and construction of a T-datum from a mutation loop.
The [API reference](api.md) records which example constructors have been
tested. Research searches and intermediate research notebooks stay in the
research checkout.

## Validation

GitHub Actions runs the regression suite, Sage doctests, an isolated wheel
installation, and all introductory notebooks in SageMath 10.8. Logs and HTML
exports are saved as workflow artifacts for inspection. The
[verification record](verification.md) distinguishes local results from
GitHub Actions results.

## Future public release

When the owner decides to make the project public:

1. Review the documented API scope and the inherited constructors listed
   as awaiting verification.
2. Confirm that the intended commit passes GitHub Actions.
3. Set the intended release version consistently in package metadata,
   `tdatum.__version__`, citation metadata, and the changelog. Update the
   wheel filename used by `tools/check_install.py` if the version changes.
4. Change repository visibility to public and update the private-access
   notices in the README and contributing guide.
5. If a versioned download is desired, tag the selected commit and attach
   its distribution files to a GitHub Release.

The current preparation stage keeps the repository private and the version
at `0.1.0.dev0`.
