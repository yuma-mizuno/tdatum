"""Build a wheel and test it in an isolated directory, without installing globally."""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def run(args, cwd, env, log):
    result = subprocess.run(args, cwd=cwd, env=env, capture_output=True, text=True, encoding="utf-8")
    log.append(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(log[-1])


def main():
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PIP_DISABLE_PIP_VERSION_CHECK="1")
    env.pop("PYTHONPATH", None)
    log = []
    dist = ROOT / "dist"
    dist.mkdir(exist_ok=True)
    run([sys.executable, "-m", "pip", "wheel", "--no-deps", "--no-build-isolation", ".", "-w", str(dist)], ROOT, env, log)
    wheel = dist / "sagemath_tdatum-0.1.0.dev0-py3-none-any.whl"
    with zipfile.ZipFile(wheel) as archive:
        assert "tdatum/t_datum.py" in archive.namelist()
        assert any(name.endswith("LICENSE") for name in archive.namelist())
    with tempfile.TemporaryDirectory(prefix="tdatum-install-") as directory:
        run([sys.executable, "-m", "pip", "install", "--no-deps", "--no-index", "--target", directory, str(wheel)], directory, env, log)
        env["PYTHONPATH"] = directory
        run([sys.executable, "-c", "import tdatum, pathlib, sys; assert pathlib.Path(tdatum.__file__).resolve().is_relative_to(pathlib.Path(sys.argv[1])); print(tdatum.__version__)", directory], directory, env, log)
        run([sys.executable, "-m", "unittest", "discover", "-s", str(ROOT / "tests"), "-v"], directory, env, log)
    output = ROOT / "test-results"
    output.mkdir(exist_ok=True)
    (output / "installation.txt").write_text("\n".join(log), encoding="utf-8")
    (output / "installation.json").write_text(json.dumps({"wheel": wheel.name, "isolated_install": "passed", "regression_tests": "passed"}, indent=2) + "\n", encoding="utf-8")
    print("Wheel build, isolated import, and installed regression tests passed.")


if __name__ == "__main__":
    main()
