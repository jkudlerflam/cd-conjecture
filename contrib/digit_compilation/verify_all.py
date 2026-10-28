#!/usr/bin/env python3
"""Rerun the new finite components; these checks are not a counterexample."""
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCRIPTS = ("parameter/checks.py", "arithmetic/checks.py",
           "integration/checks.py", "audit/gaussian_bound_checks.py")


def run(script):
    start = time.monotonic()
    result = subprocess.run([sys.executable, str(HERE / script)], cwd=ROOT,
                            env=dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1"),
                            capture_output=True, text=True)
    record = dict(script=script, passed=result.returncode == 0,
                  runtime_seconds=round(time.monotonic() - start, 3))
    if result.returncode:
        record["error"] = result.stdout + result.stderr
    return record


def main():
    start = time.monotonic()
    with ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(run, SCRIPTS))
    report = dict(scope="Local digit components and Gaussian bound; no complete core or counterexample",
                  all_passed=all(x["passed"] for x in results), tests=results,
                  runtime_seconds=round(time.monotonic() - start, 3))
    (HERE / "verification_summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not report["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
