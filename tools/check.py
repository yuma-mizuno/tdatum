"""Run regression tests and Sage doctests from a source checkout."""

import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    sage = os.environ.get("SAGE", shutil.which("sage"))
    if not sage:
        raise SystemExit("Run with sage -python, or set SAGE to the Sage executable.")
    output_dir = ROOT / "test-results"
    output_dir.mkdir(exist_ok=True)
    commands = [
        ("regression.txt", [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"]),
        ("doctests.txt", [sage, "-t", "--warn-long", "30", "--random-seed=0", "src/tdatum/t_datum.py"]),
    ]
    for filename, command in commands:
        result = subprocess.run(command, cwd=ROOT, env=env, capture_output=True,
                                text=True, encoding="utf-8")
        output = result.stdout + result.stderr
        (output_dir / filename).write_text(output, encoding="utf-8")
        print(output[-2500:])
        if result.returncode:
            raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
