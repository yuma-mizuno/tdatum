# Changelog

## 0.1.0.dev0

- Extracted the TDatum and MutationLoop code and example constructors into
  an installable SageMath package.
- Completed constructor checks for integer coefficients, a common `N0`,
  positive degrees, and a positive diagonal symmetrizer and its identities.
- Preserved symmetrizers through sign duality, mutation-loop construction,
  and loop inversion.
- Corrected loop equality to compare mutation blocks, retaining time steps.
- Validated commuting mutation blocks and loop symmetrizers; incomplete
  loops are rejected when converting to T-data.
- Preserved custom polynomial variable names in loop-to-datum conversion.
- Copied mutable inputs and protected stored matrices from mutation.
- Added regression tests, standalone Sage doctests, and introductory notebooks.
