# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Replay two dependency-free historical models in temporary copies.

This compares model outputs with saved records. It is a reproduction check, not
empirical or independent scientific validation. Historical files are not edited.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
CASES = (
    ("candidate-model", "candidate_model.py", "candidate-results.json"),
    ("diversity-persistence", "diversity_persistence.py", "diversity-persistence-results.json"),
)


def compare(expected: object, actual: object, path: str = "$*") -> int:
    """Compare all retained fields; allow roundoff in numeric calculations."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys():
            raise AssertionError(f"Object keys differ at {path}")
        return sum(compare(v, actual[k], f"{path}.{k}") for k, v in expected.items())
    if isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise AssertionError(f"Array length differs at {path}")
        return sum(compare(a, b, f"{path}[{i}]") for i, (a, b) in enumerate(zip(expected, actual)))
    if isinstance(expected, bool) or expected is None or isinstance(expected, str):
        if actual != expected or type(actual) is not type(expected):
            raise AssertionError(f"Value differs at {path}")
    elif isinstance(expected, (int, float)):
        if not isinstance(actual, (int, float)) or isinstance(actual, bool):
            raise AssertionError(f"Numeric type differs at {path}")
        if not math.isclose(expected, actual, rel_tol=1e-8, abs_tol=1e-10):
            raise AssertionError(f"Numeric value differs at {path}: {expected!r} vs {actual!r}")
    elif actual != expected:
        raise AssertionError(f"Value differs at {path}")
    return 1


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Optional new report; historical results stay unchanged.")
    args = parser.parse_args()
    checks = []
    with tempfile.TemporaryDirectory(prefix="research-replay-") as directory:
        for family, script, output in CASES:
            original = ROOT / "experiments" / family
            temporary = Path(directory) / family
            shutil.copytree(original, temporary, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
            expected = json.loads((original / output).read_text())
            process = subprocess.run(
                [sys.executable, str(temporary / script)], cwd=temporary,
                capture_output=True, text=True, check=False,
            )
            if process.returncode:
                raise RuntimeError(f"Replay failed for {family}: {process.stderr[-1200:]}")
            actual = json.loads((temporary / output).read_text())
            fields = compare(expected, actual)
            checks.append({"family": family, "matched": True, "leaf_values_compared": fields})
    report = {
        "kind": "historical_output_reproduction_check",
        "python": platform.python_version(),
        "relative_tolerance": 1e-8,
        "absolute_tolerance": 1e-10,
        "checks": checks,
        "historical_files_modified": False,
        "limit": "Tests supplied models against their saved outputs, not organisms or human benefit.",
    }
    rendered = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
