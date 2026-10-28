"""Run every new finite check with the current interpreter and record outcomes."""
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
JOBS = [
    ("circuit/verify.py", []),
    ("circuit/resource_checks.py", []),
    ("controller/checks.py", []),
    ("controller/integration_check.py", []),
    ("controller/length_counter_checks.py", []),
    ("audit/format_query_checks.py", []),
    ("audit/tiny_history_checks.py", []),
    ("local/checks.py", ["--export"]),
    ("local/port_projection_checks.py", []),
]


def run(job):
    name, args = job
    start = time.perf_counter()
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    completed = subprocess.run([sys.executable, str(HERE/name), *args],
                               cwd=ROOT, env=env, text=True, capture_output=True)
    result = dict(script=name, passed=completed.returncode == 0,
                  runtime_seconds=time.perf_counter()-start)
    if completed.returncode:
        result["error"] = completed.stdout + completed.stderr
    return result


def main():
    start = time.perf_counter()
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(run, JOBS))
    output = dict(scope="Finite component checks, not a complete counterexample",
                  all_passed=all(r["passed"] for r in results), tests=results,
                  runtime_seconds=time.perf_counter()-start)
    (HERE/"verification_summary.json").write_text(json.dumps(output, indent=2)+"\n")
    print(json.dumps(output, indent=2))
    if not output["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
