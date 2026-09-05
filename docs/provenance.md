# Source provenance

The original authorship header identifies Yuma Mizuno (2020). This package
was extracted from the local `y_system_lib` sources on 2026-09-05.
The first Git commit preserves the two input files before changes.

| Imported file | Original source | Original SHA-256 |
| --- | --- | --- |
| `src/tdatum/t_datum.py` | `y_system_lib/t_datum.py` | `126a97131735ce1df606f820ea8c132c54708e934c37bb7921b67d59e81fbb1c` |
| `src/tdatum/examples.py` | `y_system_lib/t_datum_examples.py` | `72470db17002bd0672bd1dc037968188be16385f8b0db9807cc0e2d66d73f790` |

The copies of `t_datum.py` in `Ysystem`, `continued fraction`, and the exponent
calculation directory were byte-identical to each other at extraction time
(SHA-256 `9d36e318b656efdc46034fa8067cd0a543768c3af24d3134cf272c9fd70749c9`).
The `y_system_lib` variant additionally contained `sign_dual()` and a
`check` constructor argument.

The RSG matrix pairs produced by the library and continued-fraction versions
were compared exactly for six parameter lists. They agreed for those inputs;
this is not a claim that the implementations agree for all lists.

The mathematical input-validation contract was checked against the T-datum
definition in Section 3.1 of
[Difference equations arising from cluster algebras](https://arxiv.org/abs/1912.05710v2).
The original GPL-3.0-or-later notice is preserved. The complete GPL version 3
text in `LICENSE` comes from the
[Free Software Foundation](https://www.gnu.org/licenses/gpl-3.0.txt).

The package uses Sage's own `ClusterQuiver`. It does not need the local
`my_quiver.py` copy or the older `cluster_admissible_pair` implementation.
Research searches and their generated outputs remain in the research checkout.

Packaging, documentation, tests, and the fixes listed in the changelog were
prepared with Codex assistance. This development record does not imply a
completed human review of every inherited constructor or algorithm.
