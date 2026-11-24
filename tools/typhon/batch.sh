#!/bin/bash
# The submit command supplies the recorded Slurm resources and log paths.
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
source ./environment.sh
set -u
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK" BLIS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export NUMEXPR_NUM_THREADS="$SLURM_CPUS_PER_TASK" PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export TYPHON_RESULTS="$PWD/results"
cp "$VIRTUAL_ENV/overlay.txt" logs/overlay.txt
module -t list > logs/modules.txt 2>&1

exec python - <<'PY'
import json
import os
import platform
import socket
import subprocess
import sys
from importlib.metadata import version
from pathlib import Path

manifest = json.loads(Path("manifest.json").read_text())
packages = ("numpy", "scipy", "sympy", "mpmath", "matplotlib", "threadpoolctl", "physics-tenpy", "cvxpy")
runtime = {
    "job_id": os.environ["SLURM_JOB_ID"],
    "node": socket.gethostname(), "platform": platform.platform(),
    "slurm": subprocess.check_output(["scontrol", "--version"], text=True).strip(),
    "python": sys.version, "executable": sys.executable,
    "packages": {p: version(p) for p in packages},
    "cpus_per_task": int(os.environ["SLURM_CPUS_PER_TASK"]),
}
Path("logs/runtime.json").write_text(json.dumps(runtime, indent=2) + "\n")
argv = manifest["command"]
os.chdir("src")
os.execvp(argv[0], argv)
PY
