"""Run the new finite checks; never label their success a c-d violation."""
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCRIPTS = (
    "generator/checks.py", "generator/resource_analysis.py",
    "bridge/checks.py", "bridge/integrated_checks.py",
    "bridge/local_gates_checks.py", "readout/checks.py",
    "audit/parameter_checks.py", "audit/sweep_audit.py", "core/checks.py",
)


def run(script):
    start = time.perf_counter()
    result = subprocess.run([sys.executable, str(HERE/script)], cwd=ROOT,
                            env=dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1"),
                            capture_output=True, text=True)
    record = dict(script=script, passed=result.returncode == 0,
                  runtime_seconds=round(time.perf_counter()-start, 3))
    if result.returncode:
        record["error"] = result.stdout + result.stderr
    return record


def main():
    start = time.perf_counter()
    with ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(run, SCRIPTS))
    report = dict(scope="Autonomous demonstration core and separate generator/parameter interfaces; not a high-channel counterexample",
                  all_passed=all(x["passed"] for x in results), tests=results,
                  runtime_seconds=round(time.perf_counter()-start, 3))
    (HERE/"verification_summary.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))
    if not report["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
