"""Build, execute, and export the introductory notebooks using a Sage kernel."""

import json
import os
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbconvert import HTMLExporter

ROOT = Path(__file__).resolve().parents[1]
md = nbformat.v4.new_markdown_cell
code = nbformat.v4.new_code_cell

NOTEBOOKS = {
    "01_getting_started": [
        md(r"""# Getting started with TDatum

Run this notebook with a **SageMath kernel**, after installing the package
with `sage -pip install --no-deps .` from the source checkout.

A T-datum consists of polynomial matrices $A_+,A_-$ and a positive integer
diagonal matrix $D$ satisfying the conditions in
[Mizuno, *Difference equations arising from cluster algebras*](https://arxiv.org/abs/1912.05710).
The constructor checks these conditions for the supplied matrices.
"""),
        code('''from sage.all import QQ, LaurentPolynomialRing, matrix, diagonal_matrix
from tdatum import TDatum

R = LaurentPolynomialRing(QQ, "z")
z = R.gen()
A_plus = matrix(R, [[1 + z**2, 0], [0, 1 + z**2]])
A_minus = matrix(R, [[1 + z**2, -z], [-z, 1 + z**2]])
td = TDatum(A_plus, A_minus)
td.triple()'''),
        md("The triple is $(N_0,N_+,N_-)$, where $A_+=N_0-N_+$ and $A_-=N_0-N_-$."),
        code('''N0, Np, Nm = td.triple()
assert A_plus == N0 - Np
assert A_minus == N0 - Nm
assert td.degrees() == (2, 2)
assert td.symmetrizer() == (1, 1)
{"degrees": td.degrees(), "symmetrizer": td.symmetrizer(), "permutation": td.permutation()}'''),
        md(r"""## A nonidentity symmetrizer

For the following nonsymmetric matrix $A_-$, specify $D=\operatorname{diag}(1,2)$.
Writing $A^\dagger=(A(z^{-1}))^{\mathsf T}$, the symplectic relation is
$A_+DA_-^\dagger=A_-DA_+^\dagger$.
"""),
        code('''A_minus = matrix(R, [[1 + z**2, -z], [-2*z, 1 + z**2]])
D = diagonal_matrix([1, 2])
td = TDatum(A_plus, A_minus, D)
assert A_plus * D * A_minus.transpose().subs({z: 1/z}) == A_minus * D * A_plus.transpose().subs({z: 1/z})
td.symmetrizer()'''),
        md("The sign dual exchanges $A_+$ and $A_-$ and keeps $D$. The Langlands dual conjugates the matrices by $D$ and changes the symmetrizer as described in the reference."),
        code('''dual = td.sign_dual()
assert dual.symmetrizer() == (1, 2)
assert dual.sign_dual().pair() == td.pair()
assert td.langlands_dual().langlands_dual().pair() == td.pair()
dual.pair()'''),
        md("## Validation failure\n\nThe constructor rejects a noninteger polynomial coefficient. Exact validation of a T-datum does not establish periodicity of its T/Y-system."),
        code('''try:
    TDatum(matrix(R, [[1 + z**2]]), matrix(R, [[1 + z**2 - z/2]]))
except ValueError as error:
    print(error)
else:
    raise AssertionError("A noninteger coefficient was accepted")'''),
    ],
    "02_rsg": [
        md(r"""# Reduced sine-Gordon examples

The `RSG` constructor uses an integer continued-fraction list.
Its `pair()` method constructs $(A_+,A_-)$; `t_datum()` validates the
constructed pair, with identity symmetrizer in the examples below.

This notebook verifies six explicit inputs. These computations do not prove
a statement about all parameter lists or determine a nonlinear period.
"""),
        code('''from tdatum.examples import RSG

constructor = RSG([3, 1])
td = constructor.t_datum()
assert td.degrees() == (2, 6)
td.pair()'''),
        md("The row order is the constructor's index order. Recording it makes the correspondence between matrix entries and mathematical labels explicit."),
        code('''list(zip(constructor.indices(), td.degrees()))'''),
        code('''cases = [[2,1], [3,1], [4,1], [3,2,2], [4,2,3], [3,2,2,2]]
records = []
for parameters in cases:
    datum = RSG(parameters).t_datum()
    ap, am = datum.pair()
    z = datum.variable()
    assert ap * am.transpose().subs({z: 1/z}) == am * ap.transpose().subs({z: 1/z})
    records.append((parameters, datum.size(), datum.degrees()))
records'''),
        md("A mutation loop can also be constructed from a T-datum. The following test recovers the original matrix pair from the maximal construction."),
        code('''loop = td.mutation_loop()
assert loop.is_complete()
recovered = loop.t_datum()
assert recovered.pair() == td.pair()
{"quiver_vertices": loop.size(), "time_steps": loop.length(), "mutations": loop.whole_length()}'''),
    ],
    "03_mutation_loops": [
        md(r"""# From a mutation loop to a T-datum

The displayed vertex labels are $1,2,3$. A time step mutates the commuting
vertices $1,3$, and the next time step mutates vertex $2$.

`MutationLoop` uses zero-based vertices in Python, while its Sage permutation
uses one-based values. The conversion is explicit below.
"""),
        code('''from sage.all import ZZ, matrix, Permutation
from tdatum import MutationLoop

B = matrix(ZZ, [[0,-1,0], [1,0,1], [0,-1,0]])
displayed_steps = [[1, 3], [2]]
api_steps = [[vertex - 1 for vertex in block] for block in displayed_steps]
nu = Permutation([1, 2, 3])
loop = MutationLoop(B, api_steps, nu)
assert loop.is_complete()
assert loop.b_matrix(loop.length()) == B
B'''),
        md("The two vertices in the first block commute because their exchange-matrix entry is zero. The endpoint permutation closes the exchange-matrix loop."),
        code('''assert B[0, 2] == 0
assert loop.length() == 2
assert loop.whole_length() == 3
{"displayed_steps": displayed_steps, "time_steps": loop.length(), "mutations": loop.whole_length(), "endpoint_permutation": nu}'''),
        md("The rows of the resulting T-datum are ordered by mutation points within one loop: first by time step, then by the order in that block."),
        code('''mutation_points = [(vertex, time) for time, block in enumerate(displayed_steps) for vertex in block]
td = loop.t_datum()
assert td.degrees() == (2, 2, 2)
{"row_labels": mutation_points, "degrees": td.degrees(), "symmetrizer": td.symmetrizer()}'''),
        code('''td.pair()'''),
        md("Here **time length** is 2 and **number of mutations** is 3. Closure of the exchange matrix and completeness of the loop do not by themselves prove periodicity of the nonlinear cluster variables or the T/Y-system. No such period is computed in this notebook."),
        code('''q_datum = loop.t_datum(variable_name="q")
assert str(q_datum.variable()) == "q"
q_datum.pair()'''),
    ],
}


def main():
    os.environ["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + os.environ.get("PYTHONPATH", "")
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    examples = ROOT / "examples"
    html_dir = ROOT / "docs" / "_build" / "examples"
    examples.mkdir(exist_ok=True)
    html_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    exporter = HTMLExporter()
    for name, cells in NOTEBOOKS.items():
        notebook = nbformat.v4.new_notebook(cells=cells)
        notebook.metadata.kernelspec = {"display_name": "SageMath", "language": "sage", "name": "sagemath"}
        NotebookClient(notebook, timeout=180, kernel_name="sagemath",
                       resources={"metadata": {"path": str(ROOT)}}).execute()
        nbformat.validate(notebook)
        with (examples / (name + ".ipynb")).open("w", encoding="utf-8") as output:
            nbformat.write(notebook, output)
        html, _ = exporter.from_notebook_node(notebook)
        (html_dir / (name + ".html")).write_text(html, encoding="utf-8")
        codes = [cell for cell in notebook.cells if cell.cell_type == "code"]
        errors = [output for cell in codes for output in cell.outputs if output.output_type == "error"]
        summary = {"notebook": name, "executed_cells": len(codes), "errors": len(errors)}
        summaries.append(summary)
        print(summary, flush=True)
    output = ROOT / "test-results"
    output.mkdir(exist_ok=True)
    (output / "notebooks.json").write_text(json.dumps(summaries, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
